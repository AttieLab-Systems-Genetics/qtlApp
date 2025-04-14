# Code Design Thoughts

Goal is to incorporate code in
[qtlApp/kalynn_R](https://github.com/AttieLab-Systems-Genetics/qtlApp/tree/main/kalynn_R)
into the qtlApp package as a collection of functions and modules.

This code contains considerable work on caching
as well as error checking and FST redesign.
Down the road (by the summer),
I expect to start adding new functions from other tools
I have developed, but first it helps me (really, our time)
to make the code modular.

Here is the challenge.
Naturally, code efforts have diverged,
necesitating a reset to collaborate further.
To do this well, we need to work in parallel
and on the same code base,
rather than on two different sets of code.
Goal is to use improvements in the helper functions and the app.
From a programming perspective
it makes sense to break out tasks into functions. 

The key 
[helper_functions.R](https://github.com/AttieLab-Systems-Genetics/qtlApp/blob/main/kalynn_R/modules/helper_functions.R)
have improved in cool ways.

-	Functions trait_scan, peak_finder and QTL_plot_visualizer are improvements
on Chris’s original code with caching and error checking,
as well as csv2fst and fst_rows.
-	These together are over 450 lines of code – manageable but getting long.
I prefer to see these in separate files so that they are each about
a page long (50-100 lines) and can be improved further, separately.
-	My versions are separate files in
[qltApp/R](https://github.com/AttieLab-Systems-Genetics/qtlApp/blob/main/R),
each 50-100 lines long.
I am in the process of updating them to use your improved code.

Data entry with 
[data_import.R](https://github.com/AttieLab-Systems-Genetics/qtlApp/blob/main/kalynn_R/modules/data_import.R)
is a nice organizational approach.
I have modified this (in my code) to use a “flat” CSV file
[data/import.csv](https://github.com/AttieLab-Systems-Genetics/qtlApp/blob/main/data/import.csv)
that is read in via
[importApp.R](https://github.com/AttieLab-Systems-Genetics/qtlApp/blob/main/R/importApp.R).

There is a lot of work into the server.R and ui.R,
which are of length ~1000 and ~300, respectively.
These are getting large to manage, particularly
as we both add features over time.
There are two parts of modularizing this code:

-	Pulling out R code that can be functions (notably ggplot2 code)
-	Separating out well-defined tasks into shiny modules (see below)

These extractions will simplify the server and ui considerably,
which will be important for code debugging and improvement
as we add functionality to the app.

App pieces do not come out as easily as functions,
but can be organized into shiny modules,
which I had done earlier in files in
[qltApp/R](https://github.com/AttieLab-Systems-Genetics/qtlApp/blob/main/R)
that are now called *App.R.
The goal is to have modular apps that can be separately improved
and can be interconnected like legos.
I am having to rewrite them to incorporate new code.
Here is the status:

-	importApp.R	# data import based on CSV file (working; used by all modules)
-	mainParApp.R	# parameter input via shiny (working; used by all other modules)
-	traitApp.R	# test routine to look at import files
-	scanApp.R	# wrapper for trait_scan() output (rewrite in progress)
-	peakApp.R	# wrapper for peak_finder() output (rewrite in progress)
-	qtlApp.R	# calls mainPar, scan and peak modules (more or less working)
