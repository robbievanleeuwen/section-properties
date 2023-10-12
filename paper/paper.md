---
title: "sectionproperties: A Python package for the analysis of arbitrary cross-sections using the finite element method"
tags:
  - Python
  - finite element method
  - cross-section
  - stress analysis
  - engineering
authors:
  - given-names: Robbie
    non-dropping-particle: van
    surname: Leeuwen
    orcid: 0009-0004-8056-3977
    affiliation: 1
    corresponding: true
  - name: Connor Ferster
    orcid: 0000-0000-0000-0000 # update
    affiliation: 2
affiliations:
  - name: Independent Researcher, Australia
    index: 1
  - name: Institution Name, Country
    index: 2
date: 12 October 2023
bibliography: paper.bib
---

# Summary

Summary of `sectionproperties`...

Has a clear description of the high-level functionality and purpose of the software for
a diverse, non-specialist audience been provided?

# Statement of need

Statement of need...

Does the paper have a section titled ‘Statement of need’ that clearly states what
problems the software is designed to solve, who the target audience is, and its relation
to other work?

Do the authors describe how this software compares to other commonly-used packages?

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$
\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.
$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int\_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for xxx.

For a quick reference, the following citation commands can be used:

- @pilkey -> "Author et al. (2001)"
- [@pilkey] -> "(Author et al., 2001)"
- [@pilkey; @larson] -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figures/arbitrary-section.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figures/arbitrary-section.png){ width=20% }

# Acknowledgements

We acknowledge the contributions from all the contributors to `sectionproperties`, a
list of which can be found
[here](https://github.com/robbievanleeuwen/section-properties/graphs/contributors).

# References
