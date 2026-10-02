#!/bin/bash

python -W ignore library_plot.py -params plot_parameters.csv -parts part_information.csv -designs dna_designs.csv -regulation reg_information.csv -output library.pdf
