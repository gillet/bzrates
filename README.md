# bzrates
bz-rates is a web-tool that takes into consideration the differential growth rate, b, and the plating efficiency, z, 
to calculate mutation rates from fluctuation assays. bz-rates implements two estimators: the Ma-Sandri-Sarkar maximum 
likelihood method and the generating function estimator. If you find bz-rates useful please site our paper: "bz-rates:
a web-tool to accurately estimate mutation rates from fluctuation assays.
This application was developed using Django v1.6 and Python2.7

# Running version of bz-rates
You can access a running version of bz-rates at: www.lcqb.upmc.fr/bz-rates

# Instructions
Please refer to www.lcqb.upmc.fr/bz-rates for usage instructions

# Citation
bz-rates paper is under submission process

#Add an estimator to bz-rates
You can add an estimator by modifying the view.py file and the contact.html file.

##view.py file
To ass your own estimator, 2 modification is required in the bzrates_v1/myform/view.py file.
This file contains the python functions of the estimators. 
1. You can add your own functions in it. 
2. Your custom estimator (functions) must be run in the "contact(request)" function at the indicated place (after 
the following commentary: ```#After this  you can run you own estimator```. Note that if you do so, the only web site 
will only output the restults of your own estimator. You might want to modify the if/else above or even add a choice
bar to the web site)
The results of your estimator must be saved in a list which name must be ```res```. The first element of the list must
be a string with the name of your estimator (e.g. 'myEstimator')

##contact.html file
This is the html file of the web site. Find the line 126 (```<!-- HERE ADD YOUR ESTIMATOR OUTPUT, like above-->```) and 
replace it by something like:
  ``` 
  {% elif res.0 == "myEstimator" %}
  <form name="Form1">
  <textarea readonly class="txtArea">{{ res.0 }}
  m	{{ res.1 }}
  μ	{{ res.2 }}
  </textarea>
  </form>
  {% endif %}
  ```


