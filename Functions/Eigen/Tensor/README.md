# tensor-processing

A MATLAB library for computing and processing structure tensors

# Copyright

Copyright (c) 2015 Daniel Forsberg
danne.forsberg@outlook.com

# About

This repository contains functionality for computing structure tensors
both from 2D and 3D data. Quadrature filters and monomials needed to 
do this are provided in the repository, as is functions showing how
to optimize your own quadrature filters. Note that to optimize your
own filters access to the kerngen toolbox, available at:
https://www.imt.liu.se/edu/courses/TBMI02/code/kerngen.zip
is needed. The kerngen toolbox is also useful for visualization of 
structure tensors based upon 2D data using the GOP coloring scheme.

# Setup

To use the code available in this repository, add the following 
lines to your startup.m file.

addpath('path-to-this-repository')

setup_tensor_processing_repository()

Note that this library is dependent on my matlab-utilities repository,
available from https://github.com/fordanic/matlab-utilities. If this 
repository is not available on the path then needed files can be downloaded
during setup of repository.

# Coding standard

Basic guidelines for a simple coding standard are given in the document 
"Matlab Programming Style Guidelines.pdf", available in the
matlab-utilities repository.

# Adding m-files

To create a suitable file header for new M-files, use the function 
create_new_m_function, available in the matlab-utilities repository.

# Licensing

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Source code provided in this repository is generally released under 
the GNU GENERAL PUBLIC LICENSE Version 3, if not otherwise stated.
