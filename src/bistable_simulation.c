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
// The bistable simulation engine. This engine treats   //
// the circuit in a time-independent fashion, as a      //
// system capable of 2 basis states.                    //
//                                                      //
//////////////////////////////////////////////////////////

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#ifdef GTK_GUI
  #include "callback_helpers.h"
  #include "support.h"
#else
  #define _(s) s
#endif /* def GTK_GUI */
#include "objects/QCADCell.h"
#include "simulation.h"
#include "bistable_simulation.h"
#include "custom_widgets.h"
#include "global_consts.h"

//!Options for the bistable simulation engine
//This variable is used by multiple source files
bistable_OP bistable_options = {12800, FALSE, 1e-3, 65, 12.9, 9.8e-22, 3.8e-23, 0.0, 2.0, 1000, 11.5, TRUE} ;

#ifdef GTK_GUI
extern int STOP_SIMULATION;
#else
static int STOP_SIMULATION = 0 ;
#endif /* def GTK_GUI */

// To each cell this structure is connected in order that this particular simulation engine can have its own variables. //
typedef struct {
    int number_of_neighbours;
    QCADCell **neighbours;
    int *neighbour_layer;
    double *Ek;
    double polarization;
    unsigned int clock;
    int stable;
    int cell_function;
} bistable_model;

static inline double bistable_determine_Ek (QCADCell *cell1, QCADCell *cell2, int layer_separation, bistable_OP *options);
static inline void bistable_refresh_all_Ek (int number_of_cell_layers, int *number_of_cells_in_layer, QCADCell ***sorted_cells, bistable_OP *options);

//-------------------------------------------------------------------//
// -- this is the main simulation procedure -- //
//-------------------------------------------------------------------//
simulation_data *run_bistable_simulation (int SIMULATION_TYPE, DESIGN *design, bistable_OP *options, VectorTable *pvt)
  {
  int i, j, k, q, Nix;
  int icLayers, icCellsInLayer;
  int idxMasterBitOrder = -1 ;
  int max_iterations_per_sample = ((bistable_OP *)options)->max_iterations_per_sample;
  BUS_LAYOUT_ITER bli ;
  time_t start_time, end_time;

  // optimization variables
  int number_of_cell_layers = 0, *number_of_cells_in_layer = NULL ;
  int number_of_cells;
  double clock_shift = (options->clock_high + options->clock_low)/2 + options->clock_shift;
  double clock_prefactor = (options->clock_high - options->clock_low) * options->clock_amplitude_factor;
  double four_pi_over_number_samples = 4.0 * PI / (double) options->number_of_samples;
  double two_pi_over_number_samples = 2.0 * PI / (double) options->number_of_samples;

  QCADCell ***sorted_cells = NULL ;
  bistable_model *cell_models;
  simulation_data *sim_data = NULL ;

  // simulator
  int iteration = 0;
  int stable = FALSE;
  double tolerance = ((bistable_OP *)options)->convergence_tolerance;
  double polarization_math;
  bistable_model *cell;

  double dPolarization;
  int *cell_indices;
  double clock_cache[4];

  //number of points to record in simulation results //
  //simulations can have millions of points and there is no need to plot them all //
  unsigned long int number_recorded_samples = (1 << (design->bus_layout->inputs->icUsed + 1)) * 30;
  unsigned long int record_interval;

  STOP_SIMULATION = FALSE;

  // -- get the starting time for the simulation -- //
  if((start_time = time (NULL)) < 0)
    fprintf(stderr, "Could not get start time\n");

  // Create per-layer cell arrays to be used by the engine
  simulation_inproc_data_new (design, &number_of_cell_layers, &number_of_cells_in_layer, &sorted_cells) ;

  // count the total number of cells
  number_of_cells = 0;
  for(i = 0; i < number_of_cell_layers; i++)
    number_of_cells += number_of_cells_in_layer[i];

  // allocate cell model array
  cell_models = (bistable_model*) g_malloc0( sizeof(bistable_model) * number_of_cells );

  // allocate cell indices array
    cell_indices = (int *) g_malloc0( sizeof(int) * number_of_cells );
    for( i = 0; i < number_of_cells; i++ ) {
        cell_indices[i] = i;
    }

  // build cell models
  k = 0;
  for(i = 0; i < number_of_cell_layers; i++)
    {
    for(j = 0; j < number_of_cells_in_layer[i] ; j++)
      {
        // attach the model parameters to each of the simulation cells //
        cell = cell_models + k;
        sorted_cells[i][j]->cell_model = cell;

        // copy cell function, clock, and polarization to cell model
        cell->polarization = qcad_cell_calculate_polarization(sorted_cells[i][j]);
        cell->clock = (sorted_cells[i][j]->cell_options).clock;
        cell->cell_function = (sorted_cells[i][j])->cell_function;
        cell->stable = 1; // assume stable at first; cells that aren't input or fixed are set to unstable later

        // -- Clear the model pointers so they are not dangling -- //
        cell->neighbours = NULL;
        cell->Ek = NULL;

        k++;
      }
    }

  // if we are performing a vector table simulation we consider only the activated inputs //
  if(SIMULATION_TYPE == VECTOR_TABLE)
   {
   for (Nix = 0 ; Nix < pvt->inputs->icUsed ; Nix++)
     if (!exp_array_index_1d (pvt->inputs, VT_INPUT, Nix).active_flag)
       ((bistable_model*) exp_array_index_1d (pvt->inputs, VT_INPUT, Nix).input->cell_model)->cell_function = QCAD_CELL_NORMAL ;
   }

  // write message to the command history window //
  command_history_message (_("Simulation found %d inputs %d outputs %d total cells\n"), design->bus_layout->inputs->icUsed, design->bus_layout->outputs->icUsed, number_of_cells) ;

  command_history_message(_("Starting initialization\n"));
  set_progress_bar_visible (TRUE) ;
  set_progress_bar_label (_("Bistable simulation:")) ;

    // if the number of samples is larger then the number of recorded samples then change the
    // time step to ensure only number_recorded_samples is used //
    if (number_recorded_samples >= options->number_of_samples) {
        record_interval = 1;
        number_recorded_samples = options->number_of_samples;
    } else {
        record_interval = (unsigned long int)ceil ((double)(options->number_of_samples - 1) / (double)(number_recorded_samples));
	number_recorded_samples = options->number_of_samples / record_interval + 1;
    }

  // -- Initialize the simualtion data structure -- //
  sim_data = g_malloc0 (sizeof(simulation_data));
  sim_data->number_of_traces = design->bus_layout->inputs->icUsed + design->bus_layout->outputs->icUsed;
  sim_data->number_samples = number_recorded_samples;
  sim_data->trace = g_malloc0 (sizeof (struct TRACEDATA) * sim_data->number_of_traces);

  // create and initialize the inputs into the sim data structure //
  for (i = 0; i < design->bus_layout->inputs->icUsed; i++)
    {
    sim_data->trace[i].data_labels = g_strdup (qcad_cell_get_label (QCAD_CELL (exp_array_index_1d (design->bus_layout->inputs, BUS_LAYOUT_CELL, i).cell))) ;
    sim_data->trace[i].drawtrace = TRUE;
    sim_data->trace[i].trace_function = QCAD_CELL_INPUT;
    sim_data->trace[i].data = g_malloc0 (sizeof (double) * sim_data->number_samples);
    }

  // create and initialize the outputs into the sim data structure //
  for (i = 0; i < design->bus_layout->outputs->icUsed; i++)
    {
    sim_data->trace[i + design->bus_layout->inputs->icUsed].data_labels = g_strdup (qcad_cell_get_label(QCAD_CELL(exp_array_index_1d (design->bus_layout->outputs, BUS_LAYOUT_CELL, i).cell))) ;
    sim_data->trace[i + design->bus_layout->inputs->icUsed].drawtrace = TRUE;
    sim_data->trace[i + design->bus_layout->inputs->icUsed].trace_function = QCAD_CELL_OUTPUT;
    sim_data->trace[i + design->bus_layout->inputs->icUsed].data = g_malloc0 (sizeof (double) * sim_data->number_samples);
    }

  // create and initialize the clock data //
  sim_data->clock_data = g_malloc0 (sizeof (struct TRACEDATA) * 4);

  for (i = 0; i < 4; i++)
    {
    sim_data->clock_data[i].data_labels = g_strdup_printf ("CLOCK %d", i);
    sim_data->clock_data[i].drawtrace = 1;
    sim_data->clock_data[i].trace_function = QCAD_CELL_FIXED; // Abusing the notation here

    sim_data->clock_data[i].data = g_malloc0 (sizeof (double) * sim_data->number_samples);

    if (SIMULATION_TYPE == EXHAUSTIVE_VERIFICATION)
      for (j = 0; j < sim_data->number_samples; j++)
        {
        sim_data->clock_data[i].data[j] = clock_prefactor * cos (((double)(1 << design->bus_layout->inputs->icUsed)) * (double) j * record_interval * four_pi_over_number_samples - PI * i / 2) + clock_shift ;
        sim_data->clock_data[i].data[j] = CLAMP (sim_data->clock_data[i].data[j], options->clock_low, options->clock_high) ;
        }
    else
//    if (SIMULATION_TYPE == VECTOR_TABLE)
      for (j = 0; j < sim_data->number_samples; j++)
        {
        sim_data->clock_data[i].data[j] = clock_prefactor * cos (((double)pvt->vectors->icUsed) * (double)j * record_interval * two_pi_over_number_samples - PI * i / 2) + clock_shift ;
        sim_data->clock_data[i].data[j] = CLAMP (sim_data->clock_data[i].data[j], options->clock_low, options->clock_high) ;
        }
    }

  // -- refresh all the kink energies to all the cells neighbours within the radius of effect -- //
  bistable_refresh_all_Ek (number_of_cell_layers, number_of_cells_in_layer, sorted_cells, options);

  // -- get and print the total initialization time -- //
  if((end_time = time (NULL)) < 0)
     fprintf(stderr, "Could not get end time\n");

  command_history_message("Total initialization time: %g s\n", (double)(end_time - start_time));

  command_history_message("Starting Simulation\n");

  set_progress_bar_fraction (0.0) ;

  // perform the iterations over all samples //
  for (j = 0; j < options->number_of_samples; j++)
    {
    if (j % 200 == 0)
      {
      // write the completion percentage to the command history window //
      set_progress_bar_fraction ( (float) j / (float) options->number_of_samples );
      // redraw the design if the user wants it to appear animated //
      if(options->animate_simulation)
        {
        // update the charges to reflect the polarizations so that they can be animated //
        for(icLayers = 0; icLayers < number_of_cell_layers; icLayers++)
          {
          for(icCellsInLayer = 0; icCellsInLayer < number_of_cells_in_layer[icLayers]; icCellsInLayer++)
            qcad_cell_set_polarization(sorted_cells[icLayers][icCellsInLayer],((bistable_model *)sorted_cells[icLayers][icCellsInLayer]->cell_model)->polarization);
          }
#ifdef DESIGNER
        redraw_async(NULL);
        gdk_flush () ;
#endif /* def DESIGNER */
        }
      }

    if (EXHAUSTIVE_VERIFICATION == SIMULATION_TYPE) {
        for (idxMasterBitOrder = 0, design_bus_layout_iter_first (design->bus_layout, &bli, QCAD_CELL_INPUT, &i) ; i > -1 ; design_bus_layout_iter_next (&bli, &i), idxMasterBitOrder++) {
            dPolarization = (-1 * sin (((double)(1 << idxMasterBitOrder)) * (double)j * FOUR_PI / (double)options->number_of_samples) > 0) ? 1 : -1;
            ((bistable_model *)exp_array_index_1d (design->bus_layout->inputs, BUS_LAYOUT_CELL, i).cell->cell_model)->polarization = dPolarization;
            if (0 == j % record_interval)
                sim_data->trace[i].data[j/record_interval] = dPolarization ;
        }
    } else { //    if (VECTOR_TABLE == SIMULATION_TYPE)
        for (design_bus_layout_iter_first (design->bus_layout, &bli, QCAD_CELL_INPUT, &i) ; i > -1 ; design_bus_layout_iter_next (&bli, &i)) {
            if (exp_array_index_1d (pvt->inputs, VT_INPUT, i).active_flag) {
                dPolarization = exp_array_index_2d (pvt->vectors, gboolean, (j * pvt->vectors->icUsed) / options->number_of_samples, i) ? 1 : -1;
                ((bistable_model *)exp_array_index_1d (pvt->inputs, VT_INPUT, i).input->cell_model)->polarization = dPolarization;
                if (0 == j % record_interval)
                    sim_data->trace[i].data[j/record_interval] = dPolarization ;
            }
        }
    }

    // randomize the order in which the cells are simulated to try and minimize numerical errors
    // associated with the imposed simulation order (fisher-yates algorithm)
    if( options->randomize_cells ) {
        for( i = number_of_cells; i--; ) {
            k = rand() % (i+1);

            q = cell_indices[k];
            cell_indices[k] = cell_indices[i];
            cell_indices[i] = q;
        }
    }

    // build clock cache for this iteration
    if (SIMULATION_TYPE == EXHAUSTIVE_VERIFICATION) {
        for( i = 0; i < 4; i++ ) {
            clock_cache[i] = clock_prefactor * cos (((double)(1 << design->bus_layout->inputs->icUsed)) * (double) j * four_pi_over_number_samples - PI * i / 2) + clock_shift ;
            clock_cache[i] = CLAMP (clock_cache[i], options->clock_low, options->clock_high);
            clock_cache[i] = 0.5 / clock_cache[i];
        }
    } else {
        for( i = 0; i < 4; i++ ) {
            clock_cache[i] = clock_prefactor * cos (((double)pvt->vectors->icUsed) * (double)j * two_pi_over_number_samples - PI * i / 2) + clock_shift ;
            clock_cache[i] = CLAMP (clock_cache[i], options->clock_low, options->clock_high) ;
            clock_cache[i] = 0.5 / clock_cache[i];
        }
    }

    // -- run the iteration with the given clock value -- //
    // -- iterate until the entire design has stabalized -- //
    iteration = 0;
    stable = TRUE;

    // initial iteration to determinte stability of each cell
    for( i = 0; i < number_of_cells; i++ ) {
        cell = cell_models + cell_indices[i];

        if ( cell->cell_function != QCAD_CELL_INPUT && cell->cell_function != QCAD_CELL_FIXED ) {
            polarization_math = 0;
            for (q = 0 ; q < cell->number_of_neighbours; q++)
                polarization_math += ((bistable_model*) cell->neighbours[q]->cell_model)->polarization * cell->Ek[q];

            polarization_math *= clock_cache[cell->clock];
            if( polarization_math > 1000.0 ) {
                polarization_math = 1;
            } else if( polarization_math < -1000.0 ) {
                polarization_math = -1;
            } else if( fabs( polarization_math ) > 0.001 ) {
                polarization_math /= sqrt (1 + polarization_math * polarization_math);
            }

            // If any cells polarization has changed beyond this threshold
            // then the entire circuit is assumed to have not converged.
            cell->stable = (fabs (polarization_math - cell->polarization) <= tolerance);
            stable = stable && cell->stable;

            // -- set the polarization of this cell -- //
            cell->polarization = polarization_math;
        }
    }
    // do more iterations until the entire design is stable
    while( !stable ) {
        if( iteration > max_iterations_per_sample ) {
            command_history_message( "Warning: Simulator did not fully converge, consider increasing max iterations per sample or convergence tolerance\n" );
            break;
        }

        // -- assume that the circuit is stable -- //
        stable = TRUE;
        for( i = 0; i < number_of_cells; i++ ) {
            cell = cell_models + cell_indices[i];

            if ( !cell->stable ) {
                polarization_math = 0;
                for (q = 0 ; q < cell->number_of_neighbours; q++)
                    polarization_math += ((bistable_model*) cell->neighbours[q]->cell_model)->polarization * cell->Ek[q];

                polarization_math *= clock_cache[cell->clock];
                if( polarization_math > 1000.0 ) {
                    polarization_math = 1;
                } else if( polarization_math < -1000.0 ) {
                    polarization_math = -1;
                } else if( fabs( polarization_math ) > 0.001 ) {
                    polarization_math /= sqrt (1 + polarization_math * polarization_math);
                }

                // update cell stability
                if( (fabs (polarization_math - cell->polarization) <= tolerance) ) {
                    cell->stable = 1;
                } else {
                    // If any cells polarization has changed beyond the threshold
                    // then the entire circuit is assumed to have not converged.
                    stable = FALSE;

                    // if this cell is unstable, none of its neighbours are either, so set them to unstable
                    for (q = 0 ; q < cell->number_of_neighbours; q++) {
                        if( ((bistable_model*) cell->neighbours[q]->cell_model)->cell_function != QCAD_CELL_INPUT &&
                            ((bistable_model*) cell->neighbours[q]->cell_model)->cell_function != QCAD_CELL_FIXED ) {
                            ((bistable_model*) cell->neighbours[q]->cell_model)->stable = 0;
                        }
                    }
                }

                // set the polarization of this cell
                cell->polarization = polarization_math;
            }
        }

        iteration++;
      } // while !stable

    if (0 == j % record_interval) {
        if (VECTOR_TABLE == SIMULATION_TYPE)
          for (design_bus_layout_iter_first (design->bus_layout, &bli, QCAD_CELL_INPUT, &i) ; i > -1 ; design_bus_layout_iter_next (&bli, &i))
            if (!exp_array_index_1d (pvt->inputs, VT_INPUT, i).active_flag)
              sim_data->trace[i].data[j/record_interval] = ((bistable_model *)exp_array_index_1d (pvt->inputs, VT_INPUT, i).input->cell_model)->polarization;

        // -- collect all the output data from the simulation -- //
        for (design_bus_layout_iter_first (design->bus_layout, &bli, QCAD_CELL_OUTPUT, &i) ; i > -1 ; design_bus_layout_iter_next (&bli, &i))
          sim_data->trace[design->bus_layout->inputs->icUsed + i].data[j/record_interval] = ((bistable_model *)exp_array_index_1d (design->bus_layout->outputs, BUS_LAYOUT_CELL, i).cell->cell_model)->polarization;
    }

    // -- if the user wants to stop the simulation then exit. -- //
    if(TRUE == STOP_SIMULATION)
      j = sim_data->number_samples;

    } // for number of samples

  // Free the neigbours and Ek array introduced by this simulation//
  for (i = 0; i < number_of_cell_layers; i++)
    {
    for (j = 0 ; j < number_of_cells_in_layer[i] ; j++)
      {
      g_free(((bistable_model *)sorted_cells[i][j]->cell_model)->neighbours);
      g_free(((bistable_model *)sorted_cells[i][j]->cell_model)->neighbour_layer);
      g_free(((bistable_model *)sorted_cells[i][j]->cell_model)->Ek);
      sorted_cells[i][j]->cell_model = NULL;
      }
    }

  simulation_inproc_data_free (&number_of_cell_layers, &number_of_cells_in_layer, &sorted_cells) ;

  // free cell models
  g_free( cell_models );
  g_free( cell_indices );

// -- get and print the total simulation time -- //
  if ((end_time = time (NULL)) < 0)
    fprintf(stderr, "Could not get end time\n");

  command_history_message ("Total simulation time: %g s\n", (double)(end_time - start_time));

  set_progress_bar_visible (FALSE) ;

  return sim_data;
  }//run_bistable

//-------------------------------------------------------------------//
// -- refreshes the array of Ek values for each cell in the design this is done to speed up the simulation
// since we can assume no design changes durring the simulation we can precompute all the Ek values then
// use them as necessary throughout the simulation -- //
static inline void bistable_refresh_all_Ek (int number_of_cell_layers, int *number_of_cells_in_layer, QCADCell ***sorted_cells, bistable_OP *options)
  {
  int icNeighbours = 0 ;
  bistable_model *cell_model = NULL ;
  int i, j, k, idx = 0, total_number_of_cells = 0;
  double radius_of_effect = ((bistable_OP *)options)->radius_of_effect ;

  for(i = 0; i < number_of_cell_layers; i++)
    total_number_of_cells+= number_of_cells_in_layer[i];

  // calculate the Ek for each cell //
  for(i = 0; i < number_of_cell_layers; i++)
    {
    for(j = 0 ; j < number_of_cells_in_layer[i] ; j++)
      {
      if (0 == (idx++) % 100)
        set_progress_bar_fraction((double)idx / (double)total_number_of_cells);

      cell_model = (bistable_model *)sorted_cells[i][j]->cell_model ;
      cell_model->neighbours = NULL;
      cell_model->neighbour_layer = NULL;
      cell_model->Ek = NULL;

      // select all neighbours within the provided radius //
      cell_model->number_of_neighbours = icNeighbours =
        select_cells_in_radius (sorted_cells, sorted_cells[i][j], radius_of_effect, i, number_of_cell_layers, number_of_cells_in_layer,
          ((bistable_OP *)options)->layer_separation, &(cell_model->neighbours), (int **)&(cell_model->neighbour_layer));

      if (icNeighbours > 0)
        {
        cell_model->Ek = g_malloc0 (sizeof (double) * icNeighbours);

        // ensure no memory allocation error has ocurred //
        if (NULL == cell_model->neighbours ||
            NULL == cell_model->Ek)
          exit (1);

        for (k = 0; k < icNeighbours; k++)
          //if(cell_model->neighbours[k]==NULL)printf("Null neighbour prior to passing into determine Ek for k = %d\n", k);
          // set the Ek of this cell and its neighbour //
          cell_model->Ek[k] = bistable_determine_Ek (sorted_cells[i][j], cell_model->neighbours[k], ABS (i - cell_model->neighbour_layer[k]), options);
          //printf("Ek = %e\n", cell_model->Ek[k]/1.602e-19);
        }
      }
    }
  }//refresh_all_Ek

//-------------------------------------------------------------------//
// Determines the Kink energy of one cell with respect to another this is defined as the energy of those
// cells having opposite polarization minus the energy of those two cells having the same polarization -- //
static inline double bistable_determine_Ek (QCADCell *cell1, QCADCell *cell2, int layer_separation, bistable_OP *options)
  {
  int k;
  int j;

  double distance = 0;
  double Constant = 1 / (FOUR_PI_EPSILON * options->epsilonR);
  double vertical_separation = (double)layer_separation * ((bistable_OP *)options)->layer_separation;

  static double same_polarization[4][4] =
    {{  QCHARGE_SQUAR_OVER_FOUR, -QCHARGE_SQUAR_OVER_FOUR,  QCHARGE_SQUAR_OVER_FOUR, -QCHARGE_SQUAR_OVER_FOUR },
     { -QCHARGE_SQUAR_OVER_FOUR,  QCHARGE_SQUAR_OVER_FOUR, -QCHARGE_SQUAR_OVER_FOUR,  QCHARGE_SQUAR_OVER_FOUR },
     {  QCHARGE_SQUAR_OVER_FOUR, -QCHARGE_SQUAR_OVER_FOUR,  QCHARGE_SQUAR_OVER_FOUR, -QCHARGE_SQUAR_OVER_FOUR },
     { -QCHARGE_SQUAR_OVER_FOUR,  QCHARGE_SQUAR_OVER_FOUR, -QCHARGE_SQUAR_OVER_FOUR,  QCHARGE_SQUAR_OVER_FOUR }};

  static double diff_polarization[4][4] =
    {{ -QCHARGE_SQUAR_OVER_FOUR,  QCHARGE_SQUAR_OVER_FOUR, -QCHARGE_SQUAR_OVER_FOUR,  QCHARGE_SQUAR_OVER_FOUR },
     {  QCHARGE_SQUAR_OVER_FOUR, -QCHARGE_SQUAR_OVER_FOUR,  QCHARGE_SQUAR_OVER_FOUR, -QCHARGE_SQUAR_OVER_FOUR },
     { -QCHARGE_SQUAR_OVER_FOUR,  QCHARGE_SQUAR_OVER_FOUR, -QCHARGE_SQUAR_OVER_FOUR,  QCHARGE_SQUAR_OVER_FOUR },
     {  QCHARGE_SQUAR_OVER_FOUR, -QCHARGE_SQUAR_OVER_FOUR,  QCHARGE_SQUAR_OVER_FOUR, -QCHARGE_SQUAR_OVER_FOUR }};

  double EnergyDiff = 0;
  double EnergySame = 0;

  g_assert (cell1 != NULL);
  g_assert (cell2 != NULL);
  g_assert (cell1 != cell2);

  for (k = 0; k < cell1->number_of_dots; k++)
    for (j = 0; j < cell2->number_of_dots; j++)
      {
      // determine the distance between the dots //
      distance = 1e-9 * determine_distance (cell1, cell2, k, j, vertical_separation);
      g_assert (distance != 0);

      EnergyDiff += diff_polarization[k][j] / distance;
      EnergySame += same_polarization[k][j] / distance;
      }//for other dots

 //printf("energy difference = %e energy same = %e\n", EnergyDiff/1.602e-19, EnergySame/1.602e-19);

  return Constant * (EnergyDiff - EnergySame);

  }// bistable_determine_Ek

void bistable_options_dump (bistable_OP *bistable_options, FILE *pfile)
  {
  fprintf (stderr, "bistable_options_dump:\n") ;
	fprintf (stderr, "bistable_options->number_of_samples         = %d\n",      bistable_options->number_of_samples) ;
	fprintf (stderr, "bistable_options->animate_simulation        = %s\n",      bistable_options->animate_simulation ? "TRUE" : "FALSE") ;
	fprintf (stderr, "bistable_options->convergence_tolerance     = %e\n",      bistable_options->convergence_tolerance) ;
	fprintf (stderr, "bistable_options->radius_of_effect          = %e [nm]\n", bistable_options->radius_of_effect) ;
	fprintf (stderr, "bistable_options->epsilonR                  = %e\n",      bistable_options->epsilonR) ;
	fprintf (stderr, "bistable_options->clock_high                = %e [J]\n",  bistable_options->clock_high) ;
	fprintf (stderr, "bistable_options->clock_low                 = %e [J]\n",  bistable_options->clock_low) ;
	fprintf (stderr, "bistable_options->clock_shift               = %e [J]\n",  bistable_options->clock_shift) ;
	fprintf (stderr, "bistable_options->clock_amplitude_factor    = %e\n",      bistable_options->clock_amplitude_factor) ;
	fprintf (stderr, "bistable_options->max_iterations_per_sample = %d\n",      bistable_options->max_iterations_per_sample) ;
	fprintf (stderr, "bistable_options->layer_separation          = %e [nm]\n", bistable_options->layer_separation) ;
	fprintf (stderr, "bistable_options->randomize_cells           = %s\n",      bistable_options->randomize_cells ? "TRUE" : "FALSE") ;
  }
