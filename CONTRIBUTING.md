# Contributing to CoolProp

Thank you for taking the time to contribute!

The following is a set of guidelines for contributing to CoolProp and its submodules, which are hosted at [CoolProp](https://github.com/CoolProp) on GitHub. These are mostly guidelines, not rules. Use your best judgment, and feel free to propose changes to this document in a pull request.

#### Table Of Contents

[I don't want to read this whole thing, I just have a question!!!](#i-dont-want-to-read-this-whole-thing-i-just-have-a-question)

[How Can I Contribute?](#how-can-i-contribute)
  * [Reporting Bugs](#reporting-bugs)
  * [Suggesting Enhancements](#suggesting-enhancements)
  * [Your First Code Contribution](#your-first-code-contribution)
  * [Pull Requests](#pull-requests-prs)

[Styleguides](#styleguides)
  * [Git Commit Messages](#git-commit-messages)
  * [C++ Code](https://github.com/CoolProp/CoolProp/wiki/Coding-Guidelines)

## I don't want to read this whole thing I just have a question!!!

> **Note:** Please don't file an issue to ask a question. You'll get faster results by using the resources below.

We have an official Google Group, where the community chimes in with helpful advice if you have questions, and a fairly detailed set of on-line documentation.

* [Discuss, the official CoolProp User Group](https://goo.gl/Pa7FBT)
* CoolProp documentation: [Release Version](http://www.coolprop.org) and [Development Version](http://www.coolprop.org/dev)

## How Can I Contribute?

### Reporting Bugs

This section guides you through submitting a bug report for CoolProp. Following these guidelines helps maintainers and the community understand your report, reproduce the behavior, and find related reports.

Before creating bug reports, please check [this list](#before-submitting-a-bug-report) as you might find out that you don't need to create one. When you are creating a bug report, please [include as many details as possible](#how-do-i-submit-a-good-issue-report). Fill out [the required template](ISSUE_TEMPLATE.md), the information it asks for helps us resolve issues faster.

> **Note:** If you find a **Closed** issue that seems like it is the same thing that you're experiencing, and you are using the latest CoolProp version that includes that issue resolution, open a new issue and include a link to the original issue in the body of your new one.

#### Before Submitting A Bug Report

* **Check the [FAQs file](https://github.com/CoolProp/CoolProp/blob/master/FAQ.md)** for a list of common questions and problems.
* **Check the [release](http://www.CoolProp.org) and [development](http://www.CoolProp.org) CoolProp documentation**.
* **Perform a [cursory search](https://github.com/search?q=+is%3Aissue+user%3ACoolProp)** to see if the problem has already been reported. If it has **and the issue is still open**, add a comment to the existing issue instead of opening a new one.

#### How Do I Submit A (Good) Issue Report?

Issues are tracked as [GitHub issues](https://guides.github.com/features/issues/). After you've determined the need to report an issue, create an issue on the CoolProp repository and provide the following information by filling in [the template](ISSUE_TEMPLATE.md).

Explain the problem and include additional details to help maintainers reproduce the problem:

* **Use a clear and descriptive title** for the issue to identify the problem.
* **Describe the exact steps which reproduce the problem** in as many details as possible. For example, start by explaining how you started CoolProp, e.g. which command exactly you used in the terminal, or how you started CoolProp otherwise. When listing steps, **don't just say what you did, but explain how you did it**. 
* **Provide specific examples to demonstrate the steps**. Include links to files or GitHub projects, or copy/pasteable snippets, which you use in those examples. If you're providing snippets in the issue, use [Markdown code blocks](https://help.github.com/articles/markdown-basics/#multiple-lines).
* **Describe the behavior you observed after following the steps** and point out what exactly is the problem with that behavior.
* **Explain which behavior you expected to see instead and why.**
* **Include screenshots and code segments** which show you following the described steps and clearly demonstrate the problem. 
* **If you're reporting that CoolProp crashed**, include the error message with a stack trace from the operating system.

Provide more context by answering these questions:

* **Did the problem start happening recently** (e.g. after updating to a new version of CoolProp) or was this always a problem?
* If the problem started happening recently, **can you reproduce the problem in an older version of CoolProp or on a different OS?** What's the most recent version in which the problem doesn't happen? 
* **Can you reliably reproduce the issue?** If not, provide details about how often the problem happens and under which conditions it normally happens.
* If the problem is related to working with a specific fluid, **does the problem happen for all fluids or only some?**

Include details about your configuration and environment:

* **Which version of CoolProp are you using?** You can get the exact version in Python by printing CoolProp.__version__ or by calling get_global_param_string("version") in most interfaces.  
* **What's the name and version of the OS you're using**?  
* **Which interface or CoolProp wrapper are you using (e.g. direct C++ calls, Python, Excel, Mathcad, etc.)  
* **If using Python, which version (e.g. 2.7, 3.6, etc.)**  

### Suggesting Enhancements

This section guides you through submitting an enhancement suggestion for CoolProp, including completely new features and minor improvements to existing functionality. Following these guidelines helps maintainers and the community understand your suggestion and find related suggestions.

Before creating enhancement suggestions, please check [this list](#before-submitting-an-enhancement-suggestion) as you might find out that you don't need to create one. When you are creating an enhancement suggestion, please [include as many details as possible](#how-do-i-submit-a-good-enhancement-suggestion). Fill in [the template](ISSUE_TEMPLATE.md), including the steps that you imagine you would take if the feature you're requesting existed.

#### Before Submitting An Enhancement Suggestion

* **Check the [development documentation](http://www.coolprop.org/dev)** — you might discover that the enhancement is already available or planned for the next official release. Most importantly, check if you're using [the latest version of CoolProp](http://www.coolprop.org/dev/coolprop/changelog.html) and if you can get the desired behavior by updating CoolProp.  
* **Perform a [cursory search](https://github.com/search?q=+is%3Aissue+label%3Awishlist+user%3ACoolProp)** to see if the enhancement has already been suggested. If it has, add a comment to the existing issue instead of opening a new one.

#### How Do I Submit A (Good) Enhancement Suggestion?

Enhancement suggestions are tracked as [GitHub issues](https://guides.github.com/features/issues/). Create an issue on that repository and provide the following information:

* **Use a clear and descriptive title** for the issue to identify the suggestion.
* **Provide a step-by-step description of the suggested enhancement** in as many details as possible.
* **Provide specific examples to demonstrate the steps**. Include copy/pasteable snippets which you use in those examples, as [Markdown code blocks](https://help.github.com/articles/markdown-basics/#multiple-lines).
* **Describe the current behavior** and **explain which behavior you expected to see instead** and why.
* **Include screenshots and animated GIFs** which help you demonstrate the steps or point out the part of CoolProp which the suggestion is related to. 
* **Explain why this enhancement would be useful** to most CoolProp users. 
* **Specify which version of CoolProp you're using.** 
* **Specify the name and version of the OS you're using.**

### Your First Code Contribution

Unsure where to begin contributing to CoolProp? You can start by looking through these `beginner` and `help-wanted` issues:

* [Beginner issues][beginner] - issues which should only require a few lines of code, and a test or two.
* [Help wanted issues][help-wanted] - issues which should be a bit more involved than `beginner` issues.

Both issue lists are sorted by total number of comments. While not perfect, number of comments is a reasonable proxy for impact a given change will have.

#### Local development

CoolProp can be developed locally on your machine.  Once code changes are completed and tested, make a Pull Request (PR) to the CoolProp repository.  Please see the [Wiki](https://github.com/CoolProp/CoolProp/wiki) on contributing to CoolProp.  

### Pull Requests (PRs)

* Fill in [the required template](PULL_REQUEST_TEMPLATE.md)
* Do not include issue numbers in the PR title
* Include screenshots in your pull request whenever possible.
* Follow the [C++ Coding Guidelines](https://github.com/CoolProp/CoolProp/wiki/Coding-Guidelines).
* Document new code based on the Documentation Styleguide
* Avoid platform-dependent code 

## Styleguides

### Git Commit Messages

* Use the present tense ("Add feature" not "Added feature")
* Use the imperative mood ("Move cursor to..." not "Moves cursor to...")
* Limit the first line to 72 characters or less
* Reference issues and pull requests liberally after the first line
* When only changing documentation (i.e. no actual code), include `[ci skip]` in the commit title


[beginner]:https://github.com/search?utf8=%E2%9C%93&q=is%3Aopen+is%3Aissue+label%3Abeginner+user%3Acoolprop
[help-wanted]:https://github.com/search?q=is%3Aopen+is%3Aissue+label%3Ahelp-wanted+user%3Acoolprop  
