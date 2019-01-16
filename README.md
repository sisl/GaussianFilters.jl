| Testing | Coverage | Documentation |
| :-----: | :------: | :-----------: |
| [![Build Status](https://travis-ci.org/sisl/JuliaPackageTemplate.jl.svg?branch=master)](https://travis-ci.org/sisl/JuliaPackageTemplate.jl) | [![Coverage Status](https://coveralls.io/repos/github/sisl/JuliaPackageTemplate.jl/badge.svg?branch=master)](https://coveralls.io/github/sisl/JuliaPackageTemplate.jl?branch=master) | [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://sisl.github.io/JuliaPackageTemplate.jl/latest) |

# JuliaPackageTemplate.jl
JuliaPackageTemplate provides an example Julia project template to quickly setup
continuous integration, test coverage reports, and automatic documentation deployment.

## Documentation

The documentation for the package can be found here: <https://sisl.github.io/JuliaPackageTemplate.jl/latest>

More example code and examples will be added as time permits.

## Configuring Package Name

To start off with, look through the package and replace `JuliaPackageTemplate` 
with your own project name.

`deps` contains C/C++ file dependencies of the packages which are compiled when
the package is installed by using the BinDeps.jl package.

`docs` contains

`src` contains the Julia source code of the package.

`test` contains unit tests which can be run locally

## Testing Locally

It is possible to test the package and code locally before commiting the update
and triggering a CI build. 

First, open a terminal window and navigate to the package root directory and 
start Julia

```bash
deddy@Andromeda:~$ cd /Stanford/repos/JuliaPackageTemplate.jl
deddy@Andromeda:~/Stanford/repos/JuliaPackageTemplate.jl$ julia
julia>
```

Next, activate the local package development environment 
```julia
julia> ]
(v1.0) pkg> activate .
(JuliaPackageTemplate) pkg> 
```

From here we can test the package by simply typing in the `test` command:
```julia
(JuliaPackageTemplate) pkg> test
```

If the package depends on C/C++ source files, these first must be compiled before
testing the package. In which case testing would involve two commands:
```julia
(JuliaPackageTemplate) pkg> build
(JuliaPackageTemplate) pkg> test
```

## Setting Up Continuous Integration

To setup continuous integration for the package we will use Travis-CI. Travis is 
free to use for open source projects (Thank you Travis!), or for build of private
repositories a subscription can be purchased.

To setup continuous integration your repository must contain a `.travis.yml` file 
AND continuous integration must be enabled for your repository on the TRAVIS CI 
webpage. This can be done either for your presonal repositories here:

`https://travis-ci.org/account/repositories`

or for your organizations' repositories here:

https://travis-ci.org/organizations/YOUR_GITHUB_ORGANIZATION_NAME/repositories

Note: If the project does not appear immediately, you may need to hit the "sync
repositories" button to have it appear.

## Setting Up Test Coverage

To set up test coverage, go to [coveralls.io](https://coveralls.io/repos/new),
login with your github account, and activate the project to add coverage reports.

Note: If the project does not appear immediately, you may need to hit the "sync
repositories" button to have it appear.

## Setting Up Documentation Deployment

Documentation and documentation deployment is accomplished with the Julia packager
`Documenter.jl`

To setup the automated deployment of documentation as part of the CI build process
follow the [Deploy Instructions](https://juliadocs.github.io/Documenter.jl/stable/man/hosting/).

## Adding package to Julia Package Repository

To add a package to the Julia package repository it is currently easiest to use 
[Attobot](https://github.com/attobot/attobot)
