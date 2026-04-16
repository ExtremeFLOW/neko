---
layout: default
title: High-Fidelity CFD
---

<section class="hero">
  <div class="hero-shell">
    <div class="hero-main">
      <div class="hero-brand">
        <img src="{{ '/assets/neko_cat_side_white.png' | relative_url }}" alt="Neko logo">
      </div>
      <h1>A framework for high-order fluid flow simulations</h1>
      <div class="hero-actions">
        <a class="button button-primary" href="https://github.com/ExtremeFLOW/neko/releases">Download Releases</a>
        <a class="button button-primary" href="https://github.com/ExtremeFLOW/neko">View on GitHub</a>
        <a class="button button-primary" href="https://neko.cfd/docs/release/">Documentation</a>
      </div>
    </div>
  </div>
  <div class="hero-metrics">
    <div class="metric">
      <span class="metric-value">Runs on your hardware</span>
      <ul class="metric-list">
        <li>Built with native CPUs, CUDA, HIP, and OpenCL kernels.</li>
        <li>Exceptional performance and parallel scaling.</li>
        <li>Device-abstraction layer to facilitate portable code development.</li>
      </ul>
    </div>
    <div class="metric">
      <span class="metric-value">High-order discretization</span>
      <ul class="metric-list">
        <li>Based on the spectral element method.</li>
        <li>User-selectable polynomial basis order.</li>
        <li>Hexahedral unstructured grids.</li>
      </ul>
    </div>
    <div class="metric">
      <span class="metric-value">Research ready</span>
      <ul class="metric-list">
        <li>Incompressible flow simulation with arbitrary number of additional scalars.</li>
        <li>Large-eddy simulation and wall modelling.</li>
        <li>Temporal statistics, forces, point probes, and other quality-of-life features.</li>
      </ul>
    </div>
  </div>
</section>

<section class="section">
  <div class="section-heading">
    <h2>Latest News</h2>
  </div>
  <div class="news-stack">
    {% for post in site.posts limit: 3 %}
      <a class="news-item news-item-rich" href="{{ post.url | relative_url }}">
        <span>{{ post.date | date: "%d %b %Y" }}</span>
        <strong>{{ post.title }}</strong>
        <p>{{ post.excerpt | strip_html | truncate: 150 }}</p>
      </a>
    {% endfor %}
  </div>
  <p class="section-link"><a href="{{ '/news.html' | relative_url }}">See all news</a></p>
</section>

<section class="section">
  <div class="section-heading">
    <h2>Use Cases</h2>
  </div>
  <div class="card-grid card-grid-three">
    {% for item in site.data.gallery limit: 6 %}
      <article class="media-card">
        <div class="media-card-visual media-{{ item.type }}">
          {% if item.image %}
            <img src="{{ item.image }}" alt="{{ item.title }}">
          {% elsif item.embed %}
            <iframe
              src="{{ item.embed }}"
              title="{{ item.title }}"
              loading="lazy"
              allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share"
              allowfullscreen></iframe>
          {% elsif item.url %}
            <a class="media-card-link" href="{{ item.url }}">Watch video</a>
          {% else %}
            <span>{{ item.type | upcase }}</span>
          {% endif %}
        </div>
        <div class="media-card-body">
          <h3>{{ item.title }}</h3>
          <p>{{ item.description }}</p>
        </div>
      </article>
    {% endfor %}
  </div>
  <p class="section-link"><a href="{{ '/gallery.html' | relative_url }}">Open full gallery</a></p>
</section>

<section class="section">
  <div class="section-heading">
    <h2>Publications</h2>
  </div>
  <div class="publication-list">
    {% for pub in site.data.publications limit: 4 %}
      <article class="publication-compact">
        <p class="publication-compact-year">{{ pub.year }}</p>
        <p>
          <span class="bib-authors">{{ pub.authors }}</span>
          <span class="bib-title">“{{ pub.title }}.”</span>
          <span class="bib-venue">{{ pub.venue }}.</span>
          {% if pub.doi %}<a href="{{ pub.doi }}">DOI</a>.{% endif %}
        </p>
      </article>
    {% endfor %}
  </div>
  <p class="section-link"><a href="{{ '/publications.html' | relative_url }}">Browse all publications</a></p>
</section>

<section class="section site-acknowledgements">
  <p>
    The development of Neko was supported by the European Commission Horizon 2020 project grant
    EPiGRAM-HS: Exascale Programming Models for Heterogeneous Systems (grant reference 801039), the
    European High Performance Computing Joint Undertaking (JU) and Sweden, Germany, Spain, Greece and
    Denmark under grant "CEEC - Centre of Excellence for Exascale CFD" (grant agreement No 101093393),
    the Swedish Research Council project grant Efficient Algorithms for Exascale Computational Fluid
    Dynamics (grant reference 2019-04723) and the SeRC Exascale Simulation Software Initiative (SESSI).
  </p>
  <p>
    The Neko logo was designed by Robert Hansen Jagrelius.
  </p>
  <p>
    <a class="ack-inline-brand" href="https://zulip.com/">
      <img src="https://raw.githubusercontent.com/zulip/zulip/143baa42432cde9f288bd202336ef2b11172f6e4/static/images/logo/zulip-icon-128x128.png" alt="Zulip logo">
      <span>Sponsored by Zulip</span>
    </a>, an open-source modern team chat app designed to keep both live and asynchronous conversations organized.
  </p>
</section>
