---
layout: page
title: Course Notes 
permalink: /notes/
order: 2
exclude_from_nav: false
---

<style>
.hide {
  display:none
}
</style>

<div id = 'hidden' class = 'hide' markdown="1">
* [Module 1: Course Introduction]({{ site.baseurl }}/data/notes/Intro.pdf)
* [R Script: R Basics]({{ site.baseurl }}/data/notes/R-intro.R)
* [R Script: Dataframes]({{ site.baseurl }}/data/notes/Dataframes.R)
* [Module 2: Graphical and Numerical Summaries]({{ site.baseurl }}/data/notes/Module2-Summaries.pdf)
    * [R Script: ggplot]({{ site.baseurl }}/data/notes/ggplot.R)
    * [R Script: Graphical Summaries]({{ site.baseurl }}/data/notes/graphical_summaries.R)
    * [R Script: Numerical Summaries]({{ site.baseurl }}/data/notes/numerical_summaries.R)
    * [R Script: Percentiles]({{ site.baseurl }}/data/notes/percentiles.R)
    * [Alcohol Figure]({{ site.baseurl }}/data/notes/Alcohol.pdf)
* [Module 3: Association: Contingency, Correlation, and Regression]({{ site.baseurl }}/data/notes/module3.pdf)
    * [R Script: Associations]({{ site.baseurl }}/data/notes/associations.R)
    * [R Script: Correlation and Regression]({{ site.baseurl }}/data/notes/correlation_and_regression.R)
* [Advanced R Script]({{ site.baseurl }}/data/notes/AdvancedR.R)
</div>
* [Module 4: Probability]({{ site.baseurl }}/data/notes/module4.pdf)
    * [R Script: Probability]({{ site.baseurl }}/data/notes/probability.R) 
    * [Birthday R script]({{ site.baseurl }}/data/notes/birthday.R)
    * [R Script: Normal Distribution]({{ site.baseurl }}/data/notes/normal.R) 
    * [R Script: Sampling Distributions]({{ site.baseurl }}/data/notes/sampling_distributions.R)
* [Module 5: Hypothesis Testing]({{ site.baseurl }}/data/notes/module5.pdf) 
    * [R Script: Hypothesis Testing]({{ site.baseurl }}/data/notes/hypotheses.R)  
    * [R Script: Dogs Example]({{ site.baseurl }}/data/notes/dogs.R) 
    * [R Script: Hypothesis Test for Two Sample Proportions]({{ site.baseurl }}/data/notes/two_sample_prop.R)
    * [Hypothesis Overview]({{ site.baseurl }}/data/notes/HypothesisOverview.docx) 
    * [R Script: T Distribution]({{ site.baseurl }}/data/notes/t.R) 
    * [Formula Sheet]({{ site.baseurl }}/data/notes/formula_sheet.pdf) 
    * [R Script: Two-sample t-test]({{ site.baseurl }}/data/notes/two-sample_t-test.R) 
{% comment %}
* [Module 6: Gene Expression]({{ site.baseurl }}/data/notes/GeneExpression.pdf)
    * [R Script: Microarrays]({{ site.baseurl }}/data/notes/microarray.R) 
* [Module 7: Gene Expression Analysis]({{ site.baseurl }}/data/notes/GeneExpressionAnalysis.pdf) 
    * [R Script: Contrasts]({{ site.baseurl }}/data/notes/contrasts.R)
    * [R Script: GEO and limma]({{ site.baseurl }}/data/notes/GEO-and-limma.R) 
    * [R Script: knn]({{ site.baseurl }}/data/notes/knn.R)
* [Module 8: Functional Annotation]({{ site.baseurl }}/data/notes/FunctionalAnnotation.pdf)
    * [R Script: DAVID]({{ site.baseurl }}/data/notes/DAVID.R) 
*  [Module 9: Bioinformatics Wrap-Up]({{ site.baseurl }}/data/notes/WrapUp.pdf)
    * [Movies R script]({{ site.baseurl }}/data/notes/movies.R) 
([R Notebook]({{ site.baseurl }}/data/notes/microarray.html)) 
  


{% endcomment %}
<br><br>
<center>
<div id = 'clicker'>
<a href = '#' style='font-size:120%' onclick = 'viewAll();'>Click to view all notes</a>
<script>
function viewAll() {
    document.getElementById('hidden').classList.remove('hide');
    document.getElementById('clicker').classList.add('hide');
    document.getElementsByTagName('ul')[0].style.marginBottom = '0px'
}
</script>


</div>
</center>
