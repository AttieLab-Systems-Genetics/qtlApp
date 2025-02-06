# QTL App organization

The purpose of this directory is to organize the QTL App in a modular fashion
to facilitate readability and debugging.
This is an ongoing project that is designed for QTL visualization and analysis at scale.
It is organized as a package with multiple small shiny modules, each with its own app.
The goal is to make this straightforward and easy enough for team members to develop
their own modules as the tools evolve.

To install:

```
> library(devtools)
> install_github("AttieLab-Systems-Genetics/qtlApp")
```

See
[DESCRIPTION](https://github.com/AttieLab-Systems-Genetics/qtlApp/blob/main/DESCRIPTION)
for imported packages, which are automatically installed with above step if not present.
The package uses
[roxygen2](https://roxygen2.r-lib.org/)
to build out dependencies and manual pages
([man](https://github.com/AttieLab-Systems-Genetics/qtlApp/tree/main/man))
via the command

```
> devtools::document()
```

## Callable apps

Once the package is installed, use

```
> library(qtlApp)
```

to attach the package. The following apps are available:

```
> qtlApp()     # full app (similar in appearance to previous app)
> scanApp()    # scan app (read and plot QTL genome scan)
> peakApp()    # peak app (read and show table of peaks and plot of allele estimates)
> traitApp()   # trait app (read datasets and trait file and show tables)
> mainParApp() # main parameter app (set up parameters used by other apps)
```

The app could also be deployed by

```
> shiny::runApp("inst/shinyApp/")
```

Please note that for now, the data files are hardwired in the 
[qtlSetup.R](https://github.com/AttieLab-Systems-Genetics/qtlApp/blob/main/inst/shinyApp/qtlSetup.R)
file.
Once data file location is established, this will be modified.


## Package Details

For details on modular apps, see
[ShinyApps page in Documentation repo](https://github.com/AttieLab-Systems-Genetics/Documentation/blob/main/ShinyApps.md).
The package has several analysis files used by shiny modules:

- [trait_scan.R](https://github.com/AttieLab-Systems-Genetics/qtlApp/blob/main/R/trait_scan.R)
- [QTL_plot_visualizer.R](https://github.com/AttieLab-Systems-Genetics/qtlApp/blob/main/R/QTL_plot_visualizer.R)
- [peak_finder.R](https://github.com/AttieLab-Systems-Genetics/qtlApp/blob/main/R/peak_finder.R)

The shiny modules in hierarchy of calling are:

- [qtlServer.R](https://github.com/AttieLab-Systems-Genetics/qtlApp/blob/main/R/qtlServer.R): QTL app
  - [mainParServer.R](https://github.com/AttieLab-Systems-Genetics/qtlApp/blob/main/R/mainParServer.R): main parameters
    - [traitServer.R](https://github.com/AttieLab-Systems-Genetics/qtlApp/blob/main/R/traitServer.R): break out display of `datasets` and return of `trait_list`
  - [scanServer.R](https://github.com/AttieLab-Systems-Genetics/qtlApp/blob/main/R/scanServer.R): QTL scan read from file
  - [peakServer.R](https://github.com/AttieLab-Systems-Genetics/qtlApp/blob/main/R/peakServer.R): QTL peaks read from file

The deployable app
[app.R](https://github.com/AttieLab-Systems-Genetics/qtlApp/blob/main/inst/shinyApp/app.R)
sources the file
[qtlSetup.R](https://github.com/AttieLab-Systems-Genetics/qtlApp/blob/main/inst/shinyApp/qtlSetup.R)
to load data files and then calls the 
[qtlServer.R](https://github.com/AttieLab-Systems-Genetics/qtlApp/blob/main/R/qtlServer.R).

- [qtlSetup.R](https://github.com/AttieLab-Systems-Genetics/qtlApp/blob/main/inst/shinyApp/qtlSetup.R): file setup
- [app.R](https://github.com/AttieLab-Systems-Genetics/qtlApp/blob/main/inst/shinyApp/app.R): app that calls the modules

## Plans:

The `qtlApp` is now a package on GitHub with self-documenting modules.
New modules following this design should be added as the app evolves through group discussion.

For reference on GitHub development, see

- use GitHub app on laptop (figure out tunneling to Attie Server)
- [Adding locally hosted code to GitHub](https://docs.github.com/en/migrations/importing-source-code/using-the-command-line-to-import-source-code/adding-locally-hosted-code-to-github)
- [HappyGit: Get started with GitHub](https://happygitwithr.com/usage-intro)
- [Connect your GitHub to server](https://github.com/AttieLab-Systems-Genetics/Documentation/blob/main/Server/Connecting.md#connect-your-github-to-server)
- [EDA: GitHub Pages](https://github.com/byandell-envsys/EarthDataAnalytics/blob/main/references.md#github-pages)
- [Oh Shit! Git!](https://ohshitgit.com/)

Some guiding principles:

- break out code as coherent pieces (functions or modules or workflow steps) in files.
- use indentation throughout
- use higher level packages where possible and relevant
  - be explicit about `package::function()`
  - build your own packages one at a time
- use `|>` rather than `%>%` (see 
[Differences between the base R and magrittr pipes](https://www.tidyverse.org/blog/2023/04/base-vs-magrittr-pipe/))
