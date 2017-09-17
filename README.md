# Fortran Programs

This is the collection of programs used in mechanical engineering simulations.

## Getting Started

The original source code was written in Fortran-77 dialect for Microsoft Fortran PowerStation 4.0.
At first source code was modified to be compatible with Compaq Visual Fortran 6.6.
Then more modifications allowed the programs to be run with modern gfortran-5 and
the Makefile project files were added.

Files are still encoded in [CP866](https://en.wikipedia.org/wiki/Code_page_866)
to comply with the console output in Windows systems.

## Build With

* Microsoft Fortran PowerStation 4.0
* Compaq Visual Fortran 6.6
* gfortran-5

## Structure

### Plane Deformation

Program for numerical solution of crack problem using plane deformation using
[finite element method](https://en.wikipedia.org/wiki/Finite_element_method).

### Antiplane Deformation

Program for numerical solution of crack problem with antiplane deformation using finite element method.

### Fluid

Collection of various programs calculating fluid mechanics.
