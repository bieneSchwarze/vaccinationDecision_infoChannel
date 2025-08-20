# Impact of Information Channels on the CoV-2 Vaccination Decision (part of the Study SOEP-RKI-2)
## Data Preparation and Analysis 
for a joint project between SZinn and Susanne Jordan & Sarah Jane Böttger (RKI)

For our analysis we use data from the German RKI-SOEP-2 Project ``Corona-Monitoring bundesweit" (CoMoBu-Studie)'', see https://www.rki.de/DE/Content/Gesundheitsmonitoring/Studien/lid/lid_node.html
This data is sensitive and access has to be asked for from the RKI or the SOEP. The study design is described in: Bartig, S., Brücker, H., Butschalowsky, H., Danne, C., Gößwald, A., Goßner, L., Grabka, M., . . . , Zinn, S. (2022) Corona Monitoring Nationwide (RKI-SOEP-2): Seroepidemiological Study on the Spread of SARS-CoV-2 Across Germany. Jahrbücher für Nationalökonomie und Statistik. https://doi.org/10.1515/jbnst-2022-0047.

Data access is granted on request for scientific purpose via a low-level data use request, described at the SOEP’s webpage: https://www.diw.de/en/diw_01.c.601584.en/data_access.html. Alternatively, the SOEP hotline can be contact-ed directly via soepmail@diw.de. The data access comes at no costs. 

In this Git repository you find the R source code (R version 4.1.12) used for data preparation and model estimation. 
The main file is *information_mainfile.R*. From within this file a function *fitModel()* for estimating a survey-weighted logit regression model is called. 
The respective source code is included in the file *models.R*.
Before conducting model estimation, we impute missing values by mice (package version 3.14.0) and CART. 
Estimation results are pooled by Rubin's combining rules. 
The missingness pattern in the used data set and the variables used for imputation are documented in the supplement file *SupplementMI.pdf*.
For fitting the survey-weighted logit regression models, we use Lumley's survey package (version 4.1.1) in R.
The results of a sensitivity check where the explanatory variable subjective assessment of being informed about COVID-19 has been coded differently than in the main analysis (three instead of two categories) are given in the file *sensitivityCheck_categoriesInformedness.xlsx*.
To explore the bivariate structure of the listed information channels, Figure multivariateFreq_Channels_20082025.png provides an overview. It displays the number of respondents who rated each pair among the nine information sources as very persuasive, rather persuasive, or not persuasive at all. The figure also indicates missing values.


