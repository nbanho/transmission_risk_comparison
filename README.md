# Comparing Airborne transmission risks of tuberculosis and COVID-19 in schools

This is the code accompanying the paper "Airborne transmission risks of tuberculosis and COVID-19 in schools in South Africa, Switzerland, and Tanzania: Modeling of environmental data" by Banholzer and Schmutz et al.

The repository is structured as follows:

-   *analysis*: contains code files for analyzing and visualizing the data.

-   *data-clean*: preprocessed data files

-   *data-raw*: raw data files (not yet public)

-   *preprocessing*: files to preprocess all data, i.e. environmental and epidemiological data

-   *results*: figures and tables of the paper (not in repository)

-   *utils*: additional, general code

Preprocessed (clean) data files will be made available upon publication and allow reproducing all analyses.

The relevant files regarding analysis are the following:

-   *analysis/sars.R* and *analysis/tb.R:* modeling transmission risk with a modified Wells-Riley equation, incorporating data and results from epidemiological and environmental analyses for SARS-CoV-2 and Mtb, respectively.

-   utils/quanta.R: estimation of the quanta parameter
