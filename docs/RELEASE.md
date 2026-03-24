# Release Notes Workflow (Early Stage)

easySpectra follows a simple semantic versioning cadence (`major.minor.patch`).

- Version file: root `VERSION`
- Build version: `project(... VERSION X.Y.Z)` in `CMakeLists.txt`
- Changelog: `docs/CHANGELOG.md`

## Recommended sequence

1. Implement and verify changes locally.
2. Commit in logical chunks (core/backend/gui/release).
3. Bump `VERSION` (for example `0.0.1` -> `0.0.2`).
4. Add a changelog entry with date + summary.
5. Tag the commit: `git tag v0.0.1`.
6. Push branch and tags:
   - `git push origin master`
   - `git push origin --tags`
