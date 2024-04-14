//////////////////////////////////////////////////////////
// QCADesigner                                          //
// Copyright 2002 Konrad Walus                          //
// All Rights Reserved                                  //
// Author: Konrad Walus                                 //
// Email: qcadesigner@gmail.com                         //
// WEB: http://qcadesigner.ca/                          //
//////////////////////////////////////////////////////////
//******************************************************//
//*********** PLEASE DO NOT REFORMAT THIS CODE *********//
//******************************************************//
// If your editor wraps long lines disable it or don't  //
// save the core files that way. Any independent files  //
// you generate format as you wish.                     //
//////////////////////////////////////////////////////////
// Please use complete names in variables and fucntions //
// This will reduce ramp up time for new people trying  //
// to contribute to the project.                        //
//////////////////////////////////////////////////////////
// Contents:                                            //
//                                                      //
// The coherence vector time-dependent simulation       //
// engine.                                              //
//                                                      //
//////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "objects/QCADCell.h"
#include "simulation.h"
#include "coherence_vector.h"
#include "custom_widgets.h"
#include "global_consts.h"
#ifdef GTK_GUI
  #include "callback_helpers.h"
#endif /* def GTK_GUI */

#define CONVERGENCE_TOLERANCE 1e7

// Calculates the magnitude of the 3D energy vector
//#define magnitude_energy_vector(P,G) (hypot(2*(G), (P)) * over_hbar)
#define magnitude_energy_vector(P,G) (sqrt((4.0*(G)*(G) + (P)*(P))*over_hbar_sqr)) // this way is faster, although more susceptible to overflow
// Calculates the temperature ratio
#define temp_ratio(P,G,T) (hypot((G),(P)*0.5)/((T) * kB))

#define TWO_PI 6.283185
#define FOUR_PI 12.56637061

//!Options for the coherence simulation engine
coherence_OP coherence_options = {1, 1e-15, 1e-16, 7e-11, 9.8e-22, 3.8e-23, 0.0, 2.0, 80, 12.9, 11.5, EULER_METHOD, TRUE, FALSE} ;

typedef struct {
    int number_of_neighbours;
    QCADCell **neighbours;
    int *neighbour_layer;
    double *Ek;
    double lambda_x;
    double lambda_y;
    double lambda_z;
    double polarization;
    unsigned int clock;
    int cell_function;
} coherence_model;

#ifdef GTK_GUI
extern int STOP_SIMULATION;
#else
static int STOP_SIMULATION = 0 ;
#endif /* def GTK_GUI */

// some often used variables that can be precalculated
typedef struct
  {
  double clock_prefactor;
  double clock_shift;
  double four_pi_over_number_samples;
  double two_pi_over_number_samples;
  double hbar_over_kBT;
  } coherence_optimizations;

// instance of the optimization options;
static coherence_optimizations optimization_options;

static double coherence_determine_Ek (QCADCell *cell1, QCADCell *cell2, int layer_separation, coherence_OP *options);
static void coherence_refresh_all_Ek (int number_of_cell_layers, int *number_of_cells_in_layer, QCADCell ***sorted_cells, coherence_OP *options);
static void run_coherence_iteration(int sample_number, int number_of_cells, coherence_model *cell_models, int total_number_of_inputs, unsigned long int number_samples, const coherence_OP *options, simulation_data *sim_data, int SIMULATION_TYPE, VectorTable *pvt);
static void run_coherence_stabalization(int number_of_cells, coherence_model *cell_models, int total_number_of_inputs, unsigned long int number_samples, const coherence_OP *options, simulation_data *sim_data, int SIMULATION_TYPE, VectorTable *pvt);
static inline double calculate_clock_value (unsigned int clock_num, unsigned long int sample, unsigned long int number_samples, int total_number_of_inputs, const coherence_OP *options, int SIMULATION_TYPE, VectorTable *pvt);

//-------------------------------------------------------------------//
// -- this is the main simulation procedure -- //
//-------------------------------------------------------------------//
simulation_data *run_coherence_simulation (int SIMULATION_TYPE, DESIGN *design, coherence_OP *options, VectorTable *pvt)
  {
  int i, j, k, Nix, number_of_cell_layers, *number_of_cells_in_layer;

  BUS_LAYOUT_ITER bli ;
  double dPolarization = 2.0 ;
  int idxMasterBitOrder = -1.0 ;

  int number_of_cells;
  int total_number_of_inputs = design->bus_layout->inputs->icUsed;
  time_t start_time, end_time;

  QCADCell ***sorted_cells = NULL ;
  coherence_model *cell_models;
  simulation_data *sim_data = NULL ;
  gboolean stable;

  //number of points to record in simulation results //
  //simulations can have millions of points and there is no need to plot them all //
  unsigned long int number_samples;
  unsigned long int number_recorded_samples = (1 << (design->bus_layout->inputs->icUsed + 1)) * 30;
  unsigned long int record_interval;

  STOP_SIMULATION = FALSE;

  // -- get the starting time for the simulation -- //
  if ((start_time = time (NULL)) < 0)
    fprintf (stderr, "Could not get start time\n");

  // determine the number of samples from the user options //
  number_samples = (unsigned long int)(ceil (options->duration/options->time_step));

  // if the number of samples is larger then the number of recorded samples then change the
  // time step to ensure only number_recorded_samples is used //
  if (number_recorded_samples >= number_samples)
    {
    number_recorded_samples = number_samples;
    record_interval = 1;
    }
  else
    record_interval = (unsigned long int)ceil ((double)(number_samples - 1) / (double)(number_recorded_samples));

  //fill in some of the optimizations
  optimization_options.clock_prefactor = (options->clock_high - options->clock_low) * options->clock_amplitude_factor;
  optimization_options.clock_shift = (options->clock_high + options->clock_low) * 0.5;
  optimization_options.four_pi_over_number_samples = FOUR_PI / (double)number_samples;
  optimization_options.two_pi_over_number_samples = TWO_PI / (double)number_samples;
  optimization_options.hbar_over_kBT = hbar / (kB * options->T);

  // -- spit out some messages for possible debugging -- //
  command_history_message ("About to start the coherence vector simulation with %d samples\n", number_samples);
  command_history_message ("%d samples will be recorded for graphing.\n", number_recorded_samples);
  set_progress_bar_visible (TRUE) ;
  set_progress_bar_label ("Coherence vector simulation:") ;
  set_progress_bar_fraction (0.0) ;

  // Fill in the cell arrays necessary for conducting the simulation
  simulation_inproc_data_new (design, &number_of_cell_layers, &number_of_cells_in_layer, &sorted_cells) ;

  // count the total number of cells
  number_of_cells = 0;
  for(i = 0; i < number_of_cell_layers; i++)
    number_of_cells += number_of_cells_in_layer[i];

  // allocate cell model array
  cell_models = (coherence_model*) g_malloc0( sizeof(coherence_model) * number_of_cells );

  // build cell models
  k = 0;
  for(i = 0; i < number_of_cell_layers; i++)
    for(j = 0; j < number_of_cells_in_layer[i]; j++) {
        // attach the model parameters to each of the simulation cells //
        sorted_cells[i][j]->cell_model = cell_models + k;

        // copy cell function, clock, and polarization to cell model
        ((coherence_model*) sorted_cells[i][j]->cell_model)->cell_function = sorted_cells[i][j]->cell_function;
        ((coherence_model*) sorted_cells[i][j]->cell_model)->clock = sorted_cells[i][j]->cell_options.clock;
        ((coherence_model*) sorted_cells[i][j]->cell_model)->polarization = qcad_cell_calculate_polarization( sorted_cells[i][j] );

        // -- Clear the model pointers so they are not dangling -- //
        ((coherence_model *)sorted_cells[i][j]->cell_model)->neighbours = NULL;
        ((coherence_model *)sorted_cells[i][j]->cell_model)->neighbour_layer = NULL;
        ((coherence_model *)sorted_cells[i][j]->cell_model)->Ek = NULL;

        k++;
    }

  // if we are performing a vector table simulation we consider only the activated inputs //
  if (VECTOR_TABLE == SIMULATION_TYPE)
    for (total_number_of_inputs = 0, Nix = 0 ; Nix < pvt->inputs->icUsed ; Nix++)
      {
      if (exp_array_index_1d (pvt->inputs, VT_INPUT, Nix).active_flag)
        total_number_of_inputs++ ;
      else
        // Kill the input flag for inactive inputs, so they may be correctly simulated
        ((coherence_model*) exp_array_index_1d (pvt->inputs, VT_INPUT, Nix).input->cell_model)->cell_function = QCAD_CELL_NORMAL;
      }

  // write message to the command history window //
  command_history_message ("Simulation found %d inputs %d outputs\n", total_number_of_inputs, design->bus_layout->outputs->icUsed) ;

  // build sim_data structure
  sim_data = g_malloc0 (sizeof(simulation_data)) ;
  sim_data->number_of_traces = design->bus_layout->inputs->icUsed + design->bus_layout->outputs->icUsed;
  sim_data->number_samples = number_recorded_samples;
  sim_data->trace = g_malloc0 (sizeof (struct TRACEDATA) * sim_data->number_of_traces);

  // create and initialize the inputs into the sim data structure //
  for (i = 0; i < design->bus_layout->inputs->icUsed; i++)
    {
    sim_data->trace[i].data_labels = g_strdup (qcad_cell_get_label (exp_array_index_1d (design->bus_layout->inputs, BUS_LAYOUT_CELL, i).cell));
    sim_data->trace[i].drawtrace = TRUE;
    sim_data->trace[i].trace_function = QCAD_CELL_INPUT;
    sim_data->trace[i].data = g_malloc0 (sim_data->number_samples * sizeof (double));
    }

  // create and initialize the outputs into the sim data structure //
  for (i = 0; i < design->bus_layout->outputs->icUsed; i++)
    {
    sim_data->trace[i + total_number_of_inputs].data_labels = g_strdup (qcad_cell_get_label (exp_array_index_1d (design->bus_layout->outputs, BUS_LAYOUT_CELL, i).cell));
    sim_data->trace[i + total_number_of_inputs].drawtrace = TRUE;
    sim_data->trace[i + total_number_of_inputs].trace_function = QCAD_CELL_OUTPUT;
    sim_data->trace[i + total_number_of_inputs].data = g_malloc0 (sim_data->number_samples * sizeof (double));
    }

  // create and initialize the clock data //
  sim_data->clock_data = g_malloc0 (sizeof (struct TRACEDATA) * 4);

  for (i = 0; i < 4; i++)
    {
    sim_data->clock_data[i].data_labels = g_strdup_printf ("CLOCK %d", i);
    sim_data->clock_data[i].drawtrace = 1;
    sim_data->clock_data[i].trace_function = QCAD_CELL_FIXED;
    if (NULL == (sim_data->clock_data[i].data = g_malloc0 (sim_data->number_samples * sizeof (double))))
      printf("Could not allocate memory for clock data\n");

    // fill in the clock data for the simulation results //
    for (j = 0; j<sim_data->number_samples; j++)
      //printf("j=%d, j*record_interval = %d\n",j,j*record_interval);
      sim_data->clock_data[i].data[j] = calculate_clock_value(i, j * record_interval, number_samples, total_number_of_inputs, options, SIMULATION_TYPE, pvt);
    }

    // -- refresh all the kink energies and neighbours-- //
    coherence_refresh_all_Ek (number_of_cell_layers, number_of_cells_in_layer, sorted_cells, options);

    // set initial input polarizations
    if (EXHAUSTIVE_VERIFICATION == SIMULATION_TYPE)
        for (design_bus_layout_iter_first (design->bus_layout, &bli, QCAD_CELL_INPUT, &i) ; i > -1 ; design_bus_layout_iter_next (&bli, &i))
            ((coherence_model*) exp_array_index_1d (design->bus_layout->inputs, BUS_LAYOUT_CELL, i).cell->cell_model)->polarization = sim_data->trace[i].data[0] = -1;
    else
    //  if (VECTOR_TABLE == SIMULATION_TYPE)
        for (design_bus_layout_iter_first (design->bus_layout, &bli, QCAD_CELL_INPUT, &i) ; i > -1 ; design_bus_layout_iter_next (&bli, &i))
          if (exp_array_index_1d (pvt->inputs, VT_INPUT, i).active_flag)
            ((coherence_model*) exp_array_index_1d (pvt->inputs, VT_INPUT, i).input->cell_model)->polarization = sim_data->trace[i].data[0] = exp_array_index_2d (pvt->vectors, gboolean, 0, i) ? 1 : -1;

    // Converge the steady state coherence vector for each cell so that the simulation starts without any transients //
    run_coherence_stabalization(number_of_cells, cell_models, total_number_of_inputs, number_samples, options, sim_data, SIMULATION_TYPE, pvt);

  // perform the iterations over all samples //
  for (j = 0; j < number_samples; j++)
    {
    if (0 == j % 10000)
      {
      // Update the progress bar
      set_progress_bar_fraction ((float) j / (float) number_samples) ;
      // redraw the design if the user wants it to appear animated //
#ifdef DESIGNER
      if(options->animate_simulation)
        {
        // copy all cell model polarizations to the cells
        for(i = 0; i < number_of_cell_layers; i++)
            for(k = 0; k < number_of_cells_in_layer[i]; k++)
                qcad_cell_set_polarization(sorted_cells[i][k], ((coherence_model*) sorted_cells[i][k]->cell_model)->polarization);

        redraw_async(NULL);
        gdk_flush ();
        }
#endif /* def DESIGNER */
      }

    // calculate inputs
    if (EXHAUSTIVE_VERIFICATION == SIMULATION_TYPE)
      for (idxMasterBitOrder = 0, design_bus_layout_iter_first (design->bus_layout, &bli, QCAD_CELL_INPUT, &i) ; i > -1 ; design_bus_layout_iter_next (&bli, &i), idxMasterBitOrder++)
        {
        dPolarization = (-sin (((double) (1 << idxMasterBitOrder)) * (double) j * optimization_options.four_pi_over_number_samples)) > 0 ? 1 : -1;
        ((coherence_model*) exp_array_index_1d (design->bus_layout->inputs, BUS_LAYOUT_CELL, i).cell->cell_model)->polarization = dPolarization;
        if (0 == j % record_interval)
          sim_data->trace[i].data[j/record_interval] = dPolarization ;
        }
    else
//    if (VECTOR_TABLE == SIMULATION_TYPE)
      for (design_bus_layout_iter_first (design->bus_layout, &bli, QCAD_CELL_INPUT, &i) ; i > -1 ; design_bus_layout_iter_next (&bli, &i))
        if (exp_array_index_1d (pvt->inputs, VT_INPUT, i).active_flag)
          {
            dPolarization = exp_array_index_2d (pvt->vectors, gboolean, (j*pvt->vectors->icUsed) / number_samples, i) ? 1 : -1;
            ((coherence_model*) (exp_array_index_1d( pvt->inputs, VT_INPUT, i ).input)->cell_model)->polarization = dPolarization;
          if (0 == j % record_interval)
            sim_data->trace[i].data[j/record_interval] = dPolarization ;
          }

    // run the iteration for the given sample
    run_coherence_iteration (j, number_of_cells, cell_models, total_number_of_inputs, number_samples, options, sim_data, SIMULATION_TYPE, pvt);

    // -- Set the cell model polarizations to the lambda_z value -- //
    for( k = 0; k < number_of_cells; k++ ) {
        if( (cell_models[k].cell_function == QCAD_CELL_INPUT) || (cell_models[k].cell_function == QCAD_CELL_FIXED ) )
            continue;

        if( fabs(cell_models[k].lambda_z) > 1.0 ) {
            command_history_message ("I had to abort the simulation at iteration %d because the polarization = %e was diverging.\nPossible cause is the time step is too large.\nAlternatively, you can decrease the relaxation time to reduce oscillations.\n",j, cell_models[k].lambda_z);
            command_history_message ("time step was set to %e\n", options->time_step);
            return sim_data;
        }

        cell_models[k].polarization = cell_models[k].lambda_z;
    }

    // -- collect all the output data from the simulation -- //
    if (0 == j % record_interval)
      for (design_bus_layout_iter_first (design->bus_layout, &bli, QCAD_CELL_OUTPUT, &i) ; i > -1 ; design_bus_layout_iter_next (&bli, &i))
        sim_data->trace[total_number_of_inputs + i].data[j/record_interval] = ((coherence_model*) exp_array_index_1d (design->bus_layout->outputs, BUS_LAYOUT_CELL, i).cell->cell_model)->polarization;

  if (TRUE == STOP_SIMULATION) return sim_data;

  }//for number of samples

  // Free the neigbours and Ek array introduced by this simulation//
  for (i = 0; i < number_of_cell_layers; i++)
    for (j = 0; j < number_of_cells_in_layer[i]; j++) {
      g_free (((coherence_model *)sorted_cells[i][j]->cell_model)->neighbours);
      g_free (((coherence_model *)sorted_cells[i][j]->cell_model)->neighbour_layer);
      g_free (((coherence_model *)sorted_cells[i][j]->cell_model)->Ek);
      sorted_cells[i][j]->cell_model = NULL;
    }

  simulation_inproc_data_free (&number_of_cell_layers, &number_of_cells_in_layer, &sorted_cells) ;

  // free cell models
  g_free( cell_models );

  // -- get and print the total simulation time -- //
  if ((end_time = time (NULL)) < 0)
    fprintf (stderr, "Could not get end time\n");

  command_history_message ("Total simulation time: %g s\n", (double)(end_time - start_time));
  set_progress_bar_visible (FALSE) ;
  return sim_data;
}//run_coherence


static void run_coherence_stabalization(int number_of_cells, coherence_model *cell_models, int total_number_of_inputs, unsigned long int number_samples, const coherence_OP *options, simulation_data *sim_data, int SIMULATION_TYPE, VectorTable *pvt) {
    // optimzations and options
    double options_T = options->T;

    // solver
    coherence_model *cell;
    unsigned int i, k, q;
    double clock_cache[4];
    double clock_value;
    double PEk;
    double old_lambda_x, old_lambda_y, old_lambda_z;

    int* cell_indices;
    int iterations = 0;
    int stable = FALSE;

    // allocate cell indices array
    cell_indices = (int *) g_malloc0( sizeof(int) * number_of_cells );
    for( i = 0; i < number_of_cells; i++ ) {
        cell_indices[i] = i;
    }

    // fill in clock look-up table
    for( i = 0; i < 4; i++ ) {
        clock_cache[i] = calculate_clock_value(i, 0, number_samples, total_number_of_inputs, options, SIMULATION_TYPE, pvt);
    }

    // Converge the steady state coherence vector for each cell so that the simulation starts without any transients //
    while (!stable) {
        stable = TRUE;
        iterations++;

        // randomize the order in which the cells are simulated to try and minimize numerical errors
        // associated with the imposed simulation order (fisher-yates algorithm)
        for( i = number_of_cells; i--; ) {
            k = rand() % (i+1);

            q = cell_indices[k];
            cell_indices[k] = cell_indices[i];
            cell_indices[i] = q;
        }

        for (i = 0; i < number_of_cells; i++) {
            cell = cell_models + cell_indices[i];

            if ((QCAD_CELL_INPUT == cell->cell_function) || (QCAD_CELL_FIXED == cell->cell_function))
                continue;

            PEk = 0;
            for (q = 0 ; q < cell->number_of_neighbours; q++)
                PEk += ((coherence_model*) cell->neighbours[q]->cell_model)->polarization * cell->Ek[q];

            clock_value = clock_cache[cell->clock];

            old_lambda_x = cell->lambda_x;
            old_lambda_y = cell->lambda_y;
            old_lambda_z = cell->lambda_z;

            cell->lambda_x = -2.0 * clock_value * over_hbar / magnitude_energy_vector(PEk, clock_value) * tanh(temp_ratio (PEk, clock_value, options_T));
            cell->lambda_y = 0;
            cell->lambda_z = PEk * over_hbar / magnitude_energy_vector(PEk, clock_value) * tanh(temp_ratio (PEk, clock_value, options_T));

            cell->polarization = cell->lambda_z;

            // if the lambda values are different by more then the tolerance then they have not converged //
            stable = stable && !( (fabs( cell->lambda_x - old_lambda_x ) > CONVERGENCE_TOLERANCE) ||
            (fabs( cell->lambda_y - old_lambda_y ) > CONVERGENCE_TOLERANCE) ||
            (fabs( cell->lambda_z - old_lambda_z ) > CONVERGENCE_TOLERANCE) );
        }
    }

    g_free( cell_indices );

    command_history_message ("It took %d iterations to converge the initial steady state polarization\n", iterations);
}

// -- completes one simulation iteration
// if bReportStability is flagged this also returns 1 if the lambda vectors have stabalized
static void run_coherence_iteration(int sample_number, int number_of_cells, coherence_model *cell_models, int total_number_of_inputs, unsigned long int number_samples, const coherence_OP *options, simulation_data *sim_data, int SIMULATION_TYPE, VectorTable *pvt) {
    // optimzations and options
    double optimization_options_hbar_over_kBT = optimization_options.hbar_over_kBT;
    double options_relaxation = options->relaxation;
    double options_time_step = options->time_step;

    // solver
    coherence_model *cell;
    unsigned int i, q;
    double clock_cache[4];
    double clock_value;
    double PEk;
    double lambda_x, lambda_y, lambda_z;
    double mag;
    double k1, k2, k3, k4; // for rk4
    double t1, t2; // temp

    // fill in clock look-up table
    for( i = 0; i < 4; i++ ) {
        clock_cache[i] = calculate_clock_value(i, sample_number, number_samples, total_number_of_inputs, options, SIMULATION_TYPE, pvt);
    }

    // optimization
    t1 = -options_time_step / options_relaxation;

    // loop through all the cells in the design //
    if (options->algorithm == EULER_METHOD) {
        for( i = 0; i < number_of_cells; i++ ) {
            cell = cell_models + i;

            PEk = 0;
            for (q = 0 ; q < cell->number_of_neighbours; q++)
                PEk += ((coherence_model*) cell->neighbours[q]->cell_model)->polarization * cell->Ek[q];

            // get clock for this cell
            clock_value = clock_cache[ cell->clock ];

            // Compute Magnitude
            mag = magnitude_energy_vector (PEk, clock_value);
            mag = tanh (optimization_options_hbar_over_kBT * mag) / mag; // optimization
            clock_value += clock_value;  // more optimization

            // Get lambda vector
            lambda_x = cell->lambda_x;
            lambda_y = cell->lambda_y;
            lambda_z = cell->lambda_z;

            // Calculate labmda vector
            cell->lambda_x += t1 * ((clock_value * over_hbar * mag + lambda_x) - options_relaxation * (PEk * lambda_y * over_hbar));
            cell->lambda_y += t1 * over_hbar * (options_relaxation * (PEk * lambda_x + clock_value * lambda_z) + hbar * lambda_y);
            cell->lambda_z += -t1 * over_hbar * (PEk * mag + clock_value * options_relaxation * lambda_y - hbar * lambda_z);
        }
    } else if (options->algorithm == RUNGE_KUTTA) {
         for( i = 0; i < number_of_cells; i++ ) {
            cell = cell_models + i;

            PEk = 0;
            for (q = 0 ; q < cell->number_of_neighbours; q++)
                PEk += ((coherence_model*) cell->neighbours[q]->cell_model)->polarization * cell->Ek[q];

            // Generate clock
            clock_value = clock_cache[ cell->clock ];

            // Compute magnitude
            mag = magnitude_energy_vector (PEk, clock_value);
            mag = tanh (optimization_options_hbar_over_kBT * mag) / mag; // optimization
            clock_value += clock_value;  // more optimization

            // Get lambda vector
            lambda_x = cell->lambda_x;
            lambda_y = cell->lambda_y;
            lambda_z = cell->lambda_z;

            // lambda_x
            t2 = over_hbar * (clock_value * mag - PEk * lambda_y * options_relaxation) + lambda_x;
            k1 = t1 * t2;
            k2 = t1 * (t2 + k1/2);
            k3 = t1 * (t2 + k2/2);
            k4 = t1 * (t2 + k3);
            cell->lambda_x = lambda_x + k1/6 + k2/3 + k3/3 + k4/6;

             // lambda_y
            t2 = options_relaxation * over_hbar * (PEk * lambda_x + clock_value * lambda_z) + lambda_y;
            k1 = t1 * t2;
            k2 = t1 * (t2 + k1/2);
            k3 = t1 * (t2 + k2/2);
            k4 = t1 * (t2 + k3);
            cell->lambda_y = lambda_y + k1/6 + k2/3 + k3/3 + k4/6;

            // lambda_z
            t2 = lambda_z - over_hbar * ( PEk * mag + clock_value * options_relaxation * lambda_y );
            k1 = t1 * t2;
            k2 = t1 * (t2 + k1/2);
            k3 = t1 * (t2 + k2/2);
            k4 = t1 * (t2 + k3);
            cell->lambda_z = lambda_z + k1/6 + k2/3 + k3/3 + k4/6;
        }
    } else {
        command_history_message ("coherence vector undefined algorithm\n");
        return 0;
    }
} // run_iteration


//-------------------------------------------------------------------//
// -- refreshes the array of Ek values for each cell in the design this is done to speed up the simulation
// since we can assume no design changes durring the simulation we can precompute all the Ek values then
// use them as necessary throughout the simulation -- //
static void coherence_refresh_all_Ek (int number_of_cell_layers, int *number_of_cells_in_layer, QCADCell ***sorted_cells, coherence_OP *options)
  {
  int icNeighbours = 0 ;
  coherence_model *cell_model = NULL ;
  int i,j,k;

  // calculate the Ek for each cell //
  for(i = 0 ; i < number_of_cell_layers ; i++)
    for(j = 0 ; j < number_of_cells_in_layer[i] ; j++)
      {
      cell_model = (coherence_model *)sorted_cells[i][j]->cell_model;

      // select all neighbours within the provided radius //
      cell_model->number_of_neighbours = icNeighbours =
        select_cells_in_radius(sorted_cells, sorted_cells[i][j], ((coherence_OP *)options)->radius_of_effect, i, number_of_cell_layers, number_of_cells_in_layer,
             ((coherence_OP *)options)->layer_separation, &(cell_model->neighbours), (int **)&(cell_model->neighbour_layer));

      //printf("number of neighbors = %d\n", icNeighbours);

      if (icNeighbours > 0)
        {
        cell_model->Ek = g_malloc0 (sizeof (double) * icNeighbours);

        // ensure no memory allocation error has ocurred //
        if (((coherence_model *)sorted_cells[i][j]->cell_model)->neighbours == NULL ||
            ((coherence_model *)sorted_cells[i][j]->cell_model)->Ek == NULL)
          //printf ("memory allocation error in refresh_all_Ek()\n");
          exit (1);

        for (k = 0; k < icNeighbours; k++)
          //if(cell_model->neighbours[k]==NULL)printf("Null neighbour prior to passing into determine Ek for k = %d\n", k);
          // set the Ek of this cell and its neighbour //
          cell_model->Ek[k] = coherence_determine_Ek (sorted_cells[i][j], cell_model->neighbours[k], ABS(i-cell_model->neighbour_layer[k]), options);
          //printf("Ek = %e\n", cell_model->Ek[k]/1.602e-19);
        }
      }
  }//refresh_all_Ek

//-------------------------------------------------------------------//
// Determines the Kink energy of one cell with respect to another this is defined as the energy of those
// cells having opposite polarization minus the energy of those two cells having the same polarization -- //
static double coherence_determine_Ek (QCADCell * cell1, QCADCell * cell2, int layer_separation, coherence_OP *options)
  {
  int k;
  int j;

  double distance = 0;
  double Constant = 1 / (4 * PI * EPSILON * options->epsilonR);

  double charge1[4] = { -HALF_QCHARGE,  HALF_QCHARGE, -HALF_QCHARGE,  HALF_QCHARGE };
  double charge2[4] = {  HALF_QCHARGE, -HALF_QCHARGE,  HALF_QCHARGE, -HALF_QCHARGE };

  double EnergyDiff = 0;
  double EnergySame = 0;

  g_assert (cell1 != NULL);
  g_assert (cell2 != NULL);
  g_assert (cell1 != cell2);

  for (k = 0; k < cell1->number_of_dots; k++)
    for (j = 0; j < cell2->number_of_dots; j++)
      {
      // determine the distance between the dots //
      // printf("layer seperation = %d\n", layer_seperation);
      distance = determine_distance (cell1, cell2, k, j, (double)layer_separation * ((coherence_OP *)options)->layer_separation);
      g_assert (distance != 0);

      EnergyDiff += Constant * (charge1[k] * charge2[j]) / (distance*1e-9);
      EnergySame += Constant * (charge1[k] * charge1[j]) / (distance*1e-9);
      }//for other dots

  //printf("Ek = %e\n", (EnergyDiff - EnergySame)/ 1.602e-19);

  return EnergyDiff - EnergySame;
  }// coherence_determine_Ek

//-------------------------------------------------------------------//
// Calculates the clock data at a particular sample
static inline double calculate_clock_value (unsigned int clock_num, unsigned long int sample, unsigned long int number_samples, int total_number_of_inputs, const coherence_OP *options, int SIMULATION_TYPE, VectorTable *pvt)
  {
  double clock = 0;

  if (SIMULATION_TYPE == EXHAUSTIVE_VERIFICATION)
    {
    clock = optimization_options.clock_prefactor *
      cos (((double) (1 << total_number_of_inputs)) * (double) sample * optimization_options.four_pi_over_number_samples - PI * (double)clock_num * 0.5) + optimization_options.clock_shift + options->clock_shift;

    // Saturate the clock at the clock high and low values
    clock = CLAMP (clock, options->clock_low, options->clock_high) ;
    }
  else
  if (SIMULATION_TYPE == VECTOR_TABLE)
    {
    clock = optimization_options.clock_prefactor *
      cos (((double)pvt->vectors->icUsed) * (double) sample * optimization_options.two_pi_over_number_samples - PI * (double)clock_num * 0.5) + optimization_options.clock_shift + options->clock_shift;

    // Saturate the clock at the clock high and low values
    clock = CLAMP (clock, options->clock_low, options->clock_high) ;
    }

  return clock;
  }// calculate_clock_value

void coherence_options_dump (coherence_OP *coherence_options, FILE *pfile)
  {
  fprintf (stderr, "coherence_options_dump:\n") ;
	fprintf (stderr, "coherence_options->T                         = %e [K]\n",  coherence_options->T) ;
	fprintf (stderr, "coherence_options->relaxation                = %e [s]\n",  coherence_options->relaxation) ;
	fprintf (stderr, "coherence_options->time_step                 = %e [s]\n",  coherence_options->time_step) ;
	fprintf (stderr, "coherence_options->duration                  = %e [s]\n",  coherence_options->duration) ;
	fprintf (stderr, "coherence_options->clock_high                = %e [J]\n",  coherence_options->clock_high) ;
	fprintf (stderr, "coherence_options->clock_low                 = %e [J]\n",  coherence_options->clock_low) ;
	fprintf (stderr, "coherence_options->clock_shift               = %e [J]\n",  coherence_options->clock_shift) ;
	fprintf (stderr, "coherence_options->clock_amplitude_factor    = %e\n",      coherence_options->clock_amplitude_factor) ;
	fprintf (stderr, "coherence_options->radius_of_effect          = %e [nm]\n", coherence_options->radius_of_effect) ;
	fprintf (stderr, "coherence_options->epsilonR                  = %e\n",      coherence_options->epsilonR) ;
	fprintf (stderr, "coherence_options->layer_separation          = %e [nm]\n", coherence_options->layer_separation) ;
	fprintf (stderr, "coherence_options->algorithm                 = %d\n",      coherence_options->algorithm) ;
	fprintf (stderr, "coherence_options->randomize_cells           = %s\n",      coherence_options->randomize_cells ? "TRUE" : "FALSE") ;
	fprintf (stderr, "coherence_options->animate_simulation        = %s\n",      coherence_options->animate_simulation ? "TRUE" : "FALSE") ;
  }
