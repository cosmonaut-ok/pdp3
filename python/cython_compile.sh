#!/bin/sh

cython3 --embed movie_3_component.py
gcc `python3.5-config --cflags --ldflags` -o movie_3_component movie_3_component.c
