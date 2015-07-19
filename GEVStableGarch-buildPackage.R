
# Copyrights (C) 2014 Thiago do Rego Sousa <thiagoestatistico@gmail.com>

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
#  FUNCTION:               DESCRIPTION:
#
#  Several                 Commands to create the package to send to CRAN
################################################################################



# Clean the Environment
rm(list=ls()) 

# Declare directory that contain the files
my.files.directory = '/Users/thiago/dropbox/Pappers/GEVStableGarch/GEVStableGarch/'

# Declare directory that will contain the package
cran.directory = '/Users/thiago/dropbox/Pappers/GEVStableGarch/GEVStableGarch/CRAN_versions/'

# Version of package
my.version = '1.1'

# Create a folder with this name at
# '/Users/thiago/dropbox/Pappers/GEVStableGarch/CRAN_versions/'

# Declare list of files that will be contained inside the package
# If you include/delete a new file in the package do not forget to update this list
my.list.of.files = c("class-fGEVSTABLEGARCH.R",
                     "class-fGEVSTABLEGARCHSPEC.R",
                     "dist-dskstd.R",
                     "dist-gat.R",
                     "GEVStableGarch-armaGarchDist.R",
                     "GEVStableGarch-filter.R",
                     "GEVStableGarch-fit.R",
                     "GEVStableGarch-getFormula.R",
                     "GEVStableGarch-getOrder.R",
                     "GEVStableGarch-getStart.R",
                     "GEVStableGarch-onAttach.R",
                     "GEVStableGarch-package.R",
                     "GEVStableGarch-select.R",
                     "GEVStableGarch-sim.R",
                     "GEVStableGarch-spec.R",
                     "GEVStableGarch-stationarityAparch.R",
                     "methods-show.R")

my.list.of.files.complete.adress = paste(my.files.directory, my.list.of.files, sep = '')

package.skeleton(name='GEVStableGarch',
                 path = paste(cran.directory,my.version,sep = ''),
                 code_files = my.list.of.files.complete.adress,
                 force = TRUE)

