#############################################################################################################
# Author :
#   Florian Rohart, Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created:  23-11-2015
# last modified: 18-02-2016
#
# Copyright (C) 2015
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#############################################################################################################


.onAttach <- function(libname, pkgname){ packageStartupMessage("\nLoaded mixOmics ",as.character(packageDescription("mixOmics")[["Version"]]),

    "\n\nThank you for using mixOmics! Learn how to apply our methods with our tutorials on www.mixOmics.org, vignette and bookdown on  https://github.com/mixOmicsTeam/mixOmics",
    "\nQuestions: email us at mixomics[at]math.univ-toulouse.fr  ",
    "\nBugs, Issues? https://github.com/mixOmicsTeam/mixOmics/issues",
    "\nCite us:  citation('mixOmics')"

    )}
