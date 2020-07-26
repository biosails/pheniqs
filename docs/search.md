---
layout: search
title: "Search"
permalink: /search/
id: search
---

<form action="/pheniqs/search" method="get">
  <label for="search-box">Search</label>
  <input type="text" id="search-box" name="query">
  <input type="submit" value="search">
</form>

<ul id="search-results"></ul>

<script>
  window.store = [
    {% for page in site.pages %}
        "{{ page.url | slugify }}" : {
            "id": "{{ page.url | slugify }}",
            "url": "{{ page.url | xml_escape }}",
            "title": "{{ page.title | xml_escape }}",
            "content": {{ page.content | strip_html | strip_newlines | | remove:'"' | jsonify }}
        }
      {% unless forloop.last %},{% endunless %}
    {% endfor %}
  ];
</script>
<script src="/pheniqs/js/lunr.min.js"></script>
<script src="/pheniqs/js/search.js"></script>
