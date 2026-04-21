# Neko Website

This directory contains the Jekyll-based `gh-pages` website for `neko.cfd`.
The generated API documentation published under `docs/` is separate from the Jekyll site.

## Requirements

- Ruby `3.4.9` (see `.ruby-version`)
- Bundler

## Run locally

From the repository root:

```bash
bundle config set --local path vendor/bundle
bundle install
bundle exec jekyll serve --host 127.0.0.1 --port 4000
```

Then open <http://127.0.0.1:4000>.
