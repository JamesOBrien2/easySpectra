# Release Notes Workflow (Early Stage)

For now, EasyNMR follows a very early pre-release cadence.

- Version file: root `VERSION`
- Build version: `project(... VERSION X.Y.Z)` in `CMakeLists.txt`
- Changelog: `docs/CHANGELOG.md`

## Recommended sequence

1. Implement and verify changes locally.
2. Commit in logical chunks (core/backend/gui/release).
3. Bump `VERSION` (for example `0.0.1-alpha.0` -> `0.0.1-alpha.1`).
4. Add a changelog entry with date + summary.
5. Tag the commit: `git tag v0.0.1-alpha.0`.
6. Push branch and tags:
   - `git push origin master`
   - `git push origin --tags`
