# Comparing Airborne transmission risks of tuberculosis and COVID-19 in schools

This is the code accompanying the paper "Airborne transmission risks of tuberculosis and COVID-19 in schools in South Africa, Switzerland, and Tanzania: Modeling of environmental data" by Banholzer and Schmutz et al.

The repository is structured as follows:

-   *analysis*: contains code files for analyzing and visualizing the data.

-   *data-clean*: preprocessed data files

-   *data-raw*: raw data files

-   *preprocessing*: files to preprocess all data, i.e. environmental and epidemiological data

-   *results*: figures and tables of the paper (not in repository)

-   *utils*: additional, general code

The relevant files regarding analysis are the following:

-   *analysis/infection-risk-mtb.R* and *analysis/infection-risk-sars-cov-2.R:* modeling transmission risk with a modified Wells-Riley equation, incorporating data and results from epidemiological and environmental analyses for SARS-CoV-2 and Mtb, respectively.

-   *analysis/ventilation.R*: reporting of environmental data used in the analysis.

-   utils/quanta.R: estimation of the quanta parameter.

-   utils/settings.R: Parameters in the Wells-Riley equation which have been set to fixed values and are not modeled dynamically.
