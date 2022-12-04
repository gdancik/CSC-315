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
* [Lab #5]({{ site.baseurl }}/data/hw/Lab5.R) (not collected) 
* [Lab #6]({{ site.baseurl }}/data/hw/Lab6.R) (Due: Monday, 10/24/2022)
* [Lab #7]({{ site.baseurl }}/data/hw/Lab7.R) (Due: Monday, 10/31/2022) 
* [Lab #8]({{ site.baseurl }}/data/hw/Lab8.R) (Due: Monday, 11/14/2022) 
* [Lab #9]({{ site.baseurl }}/data/hw/Lab9_contrasts.R) (Due: <span style = 'color:red'>Friday, 11/25/2022 by 5:00 PM</span><strike>Monday, 11/21/2022</strike>) 
* [Lab #10]({{ site.baseurl }}/data/hw/Lab10.R) (Due: <span style = 'color:red'>Friday, 12/02/2022</span><strike>Wednesday, 11/30/2022</strike>)
<hr>
* [Final Project]({{ site.baseurl }}/data/hw/Project.pdf) (Due: Monday, 12/12/2022 by noon; e-mail with dataset, etc due by 5:00 PM Wednesday, 12/07/2022)
    * [Example of DE genes between males/females]({{ site.baseurl }}/data/hw/SexGenes.xlsx)
    * [Real world example](https://pubmed.ncbi.nlm.nih.gov/30573692/)
{% comment %}
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
