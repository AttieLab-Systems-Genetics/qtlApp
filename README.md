# QTL App organization

The purpose of this directory is to organize the QTL App in a modular fashion
to facilitate readability and debugging.
The first step is breaking apart large files into smaller ones that are coherent.
This will soon be put on GitHub as a private repo.

Initial files:

- <appEmfinger.R>: Chris Emfinger's app emailed to Kalynn
- <appWillis.R>: Kalynn Willis's version

Setup and analysis files:

- <qtlSetup.R>: global object and library setup; source analysis functions
    - `peak_finder.R`
    - `trait_scan.R`
    - `QTL_plot_visualizer.R`

Shiny modules in hierarchy of calling:

- `qtlSetup.R`: file setup
- `app.R`: app that calls the modules
- `qtlServer.R`: QTL app components (to be broken apart later)
  - `mainParServer.R`: main parameters
    - `traitServer.R`: break out display of `datasets` and return of `trait_list`
  - `scanServer.R`: QTL scan read from file
  - `peakServer.R`: QTL peaks read from file

Plans:

- turn this folder into a GitHub repo
  - use GitHub app on laptop (figure out tunneling to Attie Server)
  - [Adding locally hosted code to GitHub](https://docs.github.com/en/migrations/importing-source-code/using-the-command-line-to-import-source-code/adding-locally-hosted-code-to-github)
  - [HappyGit: Get started with GitHub](https://happygitwithr.com/usage-intro)
  - [Connect your GitHub to server](https://github.com/AttieLab-Systems-Genetics/Documentation/blob/main/Server/Connecting.md#connect-your-github-to-server)
  - [EDA: GitHub Pages](https://github.com/byandell-envsys/EarthDataAnalytics/blob/main/references.md#github-pages)
  - [Oh Shit! Git!](https://ohshitgit.com/)
- turn repo into an R library
  - use devtools::install_github()
  - self-documenting modules

Some guiding principles:

- break out code as coherent pieces (functions or modules or workflow steps) in files.
- use indentation throughout
- use higher level packages where possible and relevant
  - be explicit about `package::function()`
  - build your own packages one at a time
- use `|>` rather than `%>%` (see 
[Differences between the base R and magrittr pipes](https://www.tidyverse.org/blog/2023/04/base-vs-magrittr-pipe/))
