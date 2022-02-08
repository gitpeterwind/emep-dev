# Research EMEP/MSC-W model(s)

## Get the source code
### I only want to run the model
Get the latest stable model version from the `master` branch
```bash
git clone -b master --single-branch git@github.com:metno/emep-mscw.git
```

#### I want a specific source/status run version
Search the [release page](https://github.com/metno/emep-mscw/releases) and note down the tag.
For example, [status 2016 version](status2016) has the tag `rv4_9`.
After cloning the repository, you can create new branch (`status2016`) from this tag as:
```bash
git checkout -b status2016 rv4_9
```

[releases]:   https://github.com/metno/emep-mscw/releases
[status2016]: https://github.com/metno/emep-mscw/releases/tag/rv4_9

### I want to develop the model
Get the develoment version from the `dev` branch
```bash
git clone -b dev --single-branch git@github.com:metno/emep-mscw.git
```

Try to:
- Commit often
- Branch your features and tasks
- Use pull requests
- Read about git [best practices][]

[cheat sheet]:    https://education.github.com/git-cheat-sheet-education.pdf
[best practices]: https://sethrobertson.github.io/GitBestPractices/

#### More git info
Print the github [cheat sheet][]

[Git Tutorial: A Comprehensive Guide][comprehensive-guide]:
- [workflow][]
- [undoing][]
- [branches][]
- [merge][] or [rebase][]

[comprehensive-guide]: https://blog.udemy.com/git-tutorial-a-comprehensive-guide/
[workflow]: https://blog.udemy.com/git-tutorial-a-comprehensive-guide/#6
[undoing]:  https://blog.udemy.com/git-tutorial-a-comprehensive-guide/#8
[branches]: https://blog.udemy.com/git-tutorial-a-comprehensive-guide/#9
[merge]:    https://blog.udemy.com/git-tutorial-a-comprehensive-guide/#10
[rebase]:   https://blog.udemy.com/git-tutorial-a-comprehensive-guide/#11
