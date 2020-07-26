---
layout: search
title: "Search"
permalink: /search/
id: search
---

<form action="/pheniqs/search.html" method="get">
  <label for="search-box">Search</label>
  <input type="text" id="search-box" name="query">
  <input type="submit" value="search">
</form>

<ul id="search-results"></ul>

<script>
  window.store = {
    {% for page in site.pages %}
        {% unless page.no_page_index %}
            "{{ page.url | slugify }}": {
              "url": "{{ page.url | xml_escape }}",
              "title": "{{ page.title | xml_escape }}",
              "content": {{ page.content | strip_html | strip_newlines | | remove:'"' | jsonify }}
            }
        {% endunless %}
      {% unless forloop.last %},{% endunless %}
    {% endfor %}
  };
</script>
<script src="/js/lunr.min.js"></script>
<script src="/js/search.js"></script>
