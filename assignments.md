---
layout: page
title: Assignments 
permalink: /assignments/
order: 3
exclude_from_nav: false
---

<style>
.due {
    background-color: yellow
}
</style>

<p style = 'color:red;font-size:104%'>Note: All assignments must be submitted through <a href = "https://easternct.blackboard.com/">Blackboard</a> unless stated otherwise. Assignments must be submitted in the correct format, which is an HTML R Notebook with none of the questions numbers modified. 
</p>
{% comment %}
{% endcomment %}

* Install <i>R/RStudio</i> and the required packages by following the instructions on the [Course Info]({{ site.baseurl }}/info/) page 
* [Lab #1]({{ site.baseurl }}/data/hw/Lab1.R) (Due: Friday, 09/09/2022) 
* [Class Survey](https://easternct.blackboard.com/) (Due: Sunday, 09/11/2022 by 5:00 PM; you may not use your grace period for this assignment)
* [Lab #2]({{ site.baseurl }}/data/hw/Lab2.R) (Due: Monday, 09/19/2022) 
* [Lab #3]({{ site.baseurl }}/data/hw/Lab3.R) (Due: Monday, 09/26/2022) 
* [Lab #4]({{ site.baseurl }}/data/hw/Lab4.R) (Due: Monday, 10/10/2022) 
<hr>
* [Lab #5]({{ site.baseurl }}/data/hw/Lab5.R) (not collected) 
{% comment %}
* [Lab #6]({{ site.baseurl }}/data/hw/Lab6.R) (Due: Monday, 10/18/2022)
* [Lab #7]({{ site.baseurl }}/data/hw/Lab7.R) (Due: Monday, 10/25/2022) 
* [Lab #8]({{ site.baseurl }}/data/hw/Lab8.R) (Due: Monday, 11/15/2022) 
* [Lab #9]({{ site.baseurl }}/data/hw/Lab9.R) (Due: Wednesday, 12/01/2022) 
* [Final Project]({{ site.baseurl }}/data/hw/Project.pdf) (Due: Wednesday, 12/08/2022 by 11:00 AM)
    * [Example of DE genes between males/females]({{ site.baseurl }}/data/hw/SexGenes.xlsx)
    * [Real world example](https://pubmed.ncbi.nlm.nih.gov/30573692/)
    * [Lab #8 results]({{ site.baseurl }}/data/hw/Lab8_results.pdf)  
* [Lab #9]({{ site.baseurl }}/data/hw/Lab9.R) (Due: Friday, 11/20/20) 
* [Classification Challenge]({{ site.baseurl }}/data/hw/Challenge.pdf) (Due: see handout)  
    * [Challenge R Script]({{ site.baseurl }}/data/hw/Challenge.R)
* Exam III survey (Due: Monday, 11/30/2020 by 11:00 AM -- see Blackboard)
    * [Lab #2 Review]({{ site.baseurl }}/data/hw/Lab2-review.R) 
    * [Lab #3 Review]({{ site.baseurl }}/data/hw/Lab3-review.R) 
    * [Lab #6 Review]({{ site.baseurl }}/data/hw/Lab6-review.R)
[[Review]({{ site.baseurl }}/data/hw/Lab6-review.R)] 
[[Solutions]({{ site.baseurl }}/data/hw/Lab7-sol.html)] 

 
{% endcomment %}

<script>
const pattern = RegExp('Due:.*([0-9]{2}/[0-9]+/[0-9]{4})');
elements = document.getElementsByTagName('li');

for (el of elements) {
        var res = pattern.exec(el.innerText);
        if (res != null && res.length >= 2) {
                if (new Date(res[1]) >= new Date()) {
                        el.className = 'due';
                }
        }
}
</script>
