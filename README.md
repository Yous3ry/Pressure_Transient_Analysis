# Pressure_Transient_Analysis
Python Based Pressure Transient Analysis (PTA) using pressure derivative plot to estimate Permeability (k), Skin (s) and wellbore storage (c)

## Assumptions:
1. Single phase (Oil) well
2. Constant production rate during flow period

## Workflow
1. Load pressure gauges data and define input parameters as a dictionary (e.g. oil rate, porosity, total compressibility, etc.).
2. Use the get_limits function to dynamically select start and end of drawdown and buildup periods (splits the data into drawdown and buildup dataframes <img align="left" width="1000" src="https://github.com/Yous3ry/Pressure_Transient_Analysis/blob/main/Pressure_Plot.png"> 
3. Use prepare data function to get the necessary parameters calculated based on test type required for analysis (drawdown or buildup)
4. Use calc_der function to calculate the second degree derivative using Bourdet derivative.
5. Use derivative_plot_analysis to draw the derivative function plot then move the horizontal line and unit slope line to estimate paratmeres. <img align="left" width="1000" src="https://github.com/Yous3ry/Pressure_Transient_Analysis/blob/main/BU_Results.png">
6. Permeability, Skin and wellbore storage are returned 
<img align="center" width="300" src="https://github.com/Yous3ry/Pressure_Transient_Analysis/blob/main/BU_Results_Numbers.png">

 
<br>
## References
1. https://www.ihsenergy.ca/support/documentation_ca/WellTest/2019_1/content/html_files/reference_materials/nomenclature.htm <br>
2. https://www.ihsenergy.ca/support/documentation_ca/WellTest/2019_1/content/html_files/analysis_types/conventional_test_analyses/derivative_analyses.htm <br>
3. https://www.ihsenergy.ca/support/documentation_ca/WellTest/content/html_files/analysis_types/conventional_test_analyses/radial_flow_analysis.htm <br>
4. https://www.ihsenergy.ca/support/documentation_ca/WellTest/2019_1/content/html_files/analysis_types/conventional_test_analyses/afterflow_analysis.htm#Summary_of_Equations_for_Afterflow_Derivative_Analysis <br>
