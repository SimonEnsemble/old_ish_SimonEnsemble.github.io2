---
layout: info
---
<center>
  School of Chemical, Biological, and Environmental Engineering<br>
</center>

<img src="osu_logo.jpg" alt="" style="width:300px">

# the professor

<a class="ppl_photo">
  <img src="photos/cory.jpg" alt="Cory Simon">
</a>

**Cory Simon** ([CV](CorySimon_academic_CV.pdf)) <br>
Assistant Professor<br>
Ph.D. Chemical Engineering. University of California, Berkeley

Hails from a small town in Ohio. Learned the ropes of scientific research at the University of Akron, Virginia Tech, Okinawa Institute of Science and Technology, University of British Columbia, Lawrence Berkeley National Laboratory, École Polytechnique Fédérale de Lausanne, and Altius Institute for Biomedical Sciences. Interned in industry at Bridgestone Research (chemical engineering) and Stitch Fix (data science).

Lives with his girlfriend, Christina, near the research forest for long walks with his Samoyed, Oslo. Digs hiking/backpacking in scenic places (photos on <a href="https://ello.co/cokes">Ello</a>), snowboarding, and wine.

Cory.Simon [at] oregonstate.edu<br>
[@CoryMSimon](https://twitter.com/CoryMSimon?lang=en)<br>
Kelley Engineering Center 2045


# the grad students

### PhD students 

{% for phd in site.data.phds %}
<a class="ppl_photo">
  <img src="{{ phd.foto }}" alt="{{ phd.name }}">
</a>
**{{ phd.name }}**. {{phd.degree }}

{{ phd.about }}

*research interests:* {{ phd.research }}

<hr>

{% endfor %}

### MS students 

{% for ms in site.data.masters %}
<a class="ppl_photo">
  <img src="{{ ms.foto }}" alt="{{ ms.name }}">
</a>
**{{ ms.name }}**. {{ms.degree }}

{{ ms.about }}

*research interests:* {{ ms.research }}

<hr>

{% endfor %}

# the undergrad students

{% for ugrad in site.data.ugrads %}
<a class="ppl_photo">
  <img src="{{ ugrad.foto }}" alt="{{ ugrad.name }}">
</a>
**{{ ugrad.name }}**. *{{ ugrad.major }}*.
{{ ugrad.about }}

*research interests:* {{ ugrad.research }}

<hr>

{% endfor %}

# alumni

## graduate students

{% for alum in site.data.grad_alumni %}
**{{ alum.name }}**. *{{ alum.major }}*.<br>
thesis: "{{ alum.thesis }}"<br>
{% if alum.gs != "missing" %}<a href="{{ alum.gs }}">Google Scholar profile</a>.
{% endif %}

{% endfor %}

## undergraduate students

{% for alum in site.data.ugrad_alumni %}
**{{ alum.name }}**. major: *{{ alum.major }}*.
{% if alum.now != "missing" %} <br> now: {{ alum.now }}
{% endif %}
{% endfor %}
