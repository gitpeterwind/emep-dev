# Branches

At stated on the [Readme][], the emep-mscw repository has 2 main branches dedicated to the research version of the EMEP/MSC-W CTM:

- *master*: Release branch. Meant for code that just "works out of the box".
- *dev*: Development branch. Bleeding edge, might be broken.

Other branches:

- *ecosx*: box model built on top of *esx*.
- *dev-3DVar*: 3DVar development.

## Update *master* branch

In order to keep the *master* stable with a linear history, we ask for our developer team to:
- Update the *master* branch via [pull request][ghPR] from the *dev* branch or a feature specific branch.
- [Merge][] the pull request using the [squash][] option

### Checklist

- Documentation and [user guide](http://emep-ctm.readthedocs.io/en/latest/)
  - Update our user guide with the changes that affect users [here](https://github.com/metno/emep-ctm/blob/docs/Input.rst) (or in the other rst files). The web pages will update automatically.
- Driver scripts (eg *run.pl*)
  - Keep the defaults unchanged, whenever possible.
  - Make sure new files/data paths are available on *vilje* and *stallo*
- Source code
  - Compile and test with full debug options (make DEBUG=yes)
  - Re-compile normally, and run a benchmark and submit to [AeroCom][]
  - Describe the update in [Log.changes][]

## Examples

```bash
# Clone only the *dev* branch
git clone -b dev --single-branch git@github.com:metno/emep-mscw.git

# add ecosx (remote) branch to (local) dev repo
cd emep-mscw/
git remote set-branches --add origin ecosx
git fetch
git checkout ecosx

# status of all branches
git branch -a -vv
```

```
* dev                  f7fd1af [origin/dev] extending .gitignore with eclipse-startup files
  ecosx                1838931 [origin/ecosx] ESX removed old DIFFS file
  remotes/origin/dev   f7fd1af extending .gitignore with eclipse-startup files
  remotes/origin/ecosx 1838931 ESX removed old DIFFS file
```

[Readme]: https://github.com/metno/emep-mscw/blob/master/README.md
[Log.changes]: https://github.com/metno/emep-mscw/blob/master/Log.changes
[AeroCom]:  http://aerocom.met.no/cgi-bin/aerocom/surfobs_annualrs.pl?PROJECT=EMEP&MODELLIST=EMEP
[ghPR]:   https://help.github.com/articles/about-pull-requests/
[merge]:  https://help.github.com/articles/merging-a-pull-request/
[squash]: https://help.github.com/articles/about-pull-request-merges/#squash-and-merge-your-pull-request-commits

# Benckmarks
Following the naming convention and procedure described on #9 and #37
- `rv4_12` (5fceaeb): vilje:~alvarov/work/Benchmark/EECCA.2012/20170420_dev@5fceaeb/
- `rv4_13` (10d9e3f): vilje:~alvarov/work/Benchmark/EECCA.2012/20170426_dev@10d9e3f/
- `rv4_15` (924fe08): vilje:~alvarov/work/Benchmark/EECCA.2012/20170606_dev@924fe08/
- `rv4_15_1` (494ad81, status 2017): vilje:~alvarov/work/Benchmark/EECCA.2012/20170828_dev@494ad81/
