# Contributing

Thank you for your interest in contributing to Mantel! This project is in the stable/maintenance phase and is not under active development. I am only aiming to fix bugs and keep it compatible with current versions of Python so that it remains useful to the community. With that in mind, contributions are very welcome, so long as they align with the guidelines below.


## Scope and Expectations

* **Maintenance-mode:** The package is not being actively expanded. I am happy to consider bug fixes, small improvements, and updates for compatibility with new Python releases.
* **No major redesigns:** Please avoid proposing large structural refactors or redesigns based on personal preference.
* **Avoid breaking changes:** The package is used in research contexts, so please do not submit changes that break the existing public API or alter core behavior in ways that may affect other researchers.

If you're unsure whether a change is appropriate, opening an issue first is strongly encouraged. This helps avoid lost effort on pull requests that cannot be merged due to scope or compatibility requirements.


## Bug Reports & Issues

If you encounter a bug, please open an issue with:

* A short, clear description of the problem.
* A minimal working example that shows the expected vs. actual behavior.
* If possible, a proposal about how to fix the problem.

If the bug is particularly straight-forward and uncontroversial, you can skip directly to the pull request.


## Pull Requests

Before starting a PR, consider whether it might be more appropriate to open an issue first to confirm that the contribution fits within the project's scope.

When submitting a PR, please include:

* A clear description of **what the PR does** and **why**, linking to the relevant issue (if applicable).
* If fixing a bug, provide a minimal example so that I can quickly understand and verify the solution.
* Ideally, please add new test code to cover new behaviors and avoid future regressions.

As always, please try to follow good "git hygiene":

* One logical change per commit.
* Don't bundle unrelated changes into a single commit or PR.
* Use concise, descriptive commit messages.
* Each commit should leave the package in a working state that passes the unit tests.

Finally, before submitting a PR, please ensure that:

1. The unit tests still pass by running `pytest tests`.
2. The code is formatted in [Black](https://black.readthedocs.io) style: `black mantel`.
3. New code follows existing conventions in the codebase.


## Thank You

Your contributions – whether through bug reports, documentation improvements, or pull requests – help keep this software useful to the research community.
