# antigenic-maps
Code to infer antigenic maps from titer panel data.

This repo contains code to infer antigenic maps from titers panel data, building on the methods of [Smith et al. 2004](), and [Bedford et al. 2014]().

The overall goal of this project was to simulate multi-epitope antigens and explore how maps fitted to synthetic titers estimated for these hypothetical multi-epitope antigens compared to epitope-specific maps.

Eventually, this project morphed into the ferret-human distance comparison, but this code could still be used for antigenic cartography.

# directories

* R - scripts containing useful functions and code.
* experiments - R notebooks designed to test code and run simulation studies. Notes within each notebook explain the purpose and approach.
* stan - .stan models used to estimate antigenic maps using HMC.
* tests-and-development - R notebooks used to test the inference code and make improvements.

# resources

* This [Asana project](https://app.asana.com/0/1201886272028192/board) contains some relevant tasks and notes.
* This [lab meeting presentation](https://app.asana.com/0/1201886272028192/1202864971519590) and [check-in slides](https://app.asana.com/0/1201886272028192/1202953282623929) outline the approach and preliminary results.
* This [writeup](https://app.asana.com/0/1201886272028192/1202042948298325) describes the model I was using to generate synthetic titers to multi-epitope antigens.
