# Contributing to PSOPT

Thank you for your interest in contributing to PSOPT! This guide outlines how to contribute code, report bugs, improve documentation, and help maintain the integrity of the project using our automated CI system.

---

## 🚀 Getting Started

1. **Fork the repository**  
   Click the “Fork” button on GitHub and clone your fork locally:
   ```bash
   git clone https://github.com/your-username/psopt.git
   cd psopt
   ```

2. **Create a branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

3. **Install dependencies**  
   Use your preferred Linux distribution, or run PSOPT inside Docker.  
   Refer to the `Dockerfile` in the repository for setup.

---

## 🧪 Running Tests Locally

Each example in `./examples/` produces an output file:  
```
psopt_solution_<example>.txt
```

Before submitting a pull request:
- Run selected examples (e.g., `twoburn`, `shuttle_reentry`, etc.)
- Confirm that the final solution file includes:
  - `NLP solver reports: The problem has been solved!`
  - `Optimal (unscaled) cost function value: ...`

Compare cost values to known reference results.

---

## 🤖 Continuous Integration (CI)

All changes are automatically validated through GitHub Actions.

Our CI performs:
- ✅ Build of the PSOPT library using Docker (Arch Linux)
- ✅ Execution of selected example problems
- ✅ Parsing of output to verify optimal cost
- ✅ Badge and dashboard generation (on `master` only)

See:
- [CI Badge](https://psopt.github.io/psopt/artifacts/examples_badge.json)
- [Dashboard](https://psopt.github.io/psopt/artifacts/index.html)

Pull requests are **not allowed to update GitHub Pages artifacts**—only merged commits to `master` will.

---

## 📌 CI Requirements for Pull Requests

- ✅ Must be based on the latest `master`
- ✅ Must **pass the Build-and-test** GitHub Action
- ❌ CI will fail if solution files are missing or invalid
- ❌ Artifacts are **not published** for PRs

---

## ✅ Commit & PR Guidelines

- Use clear, concise commit messages (e.g., `fix: handle edge case in low_thrust`)
- Separate unrelated changes into different commits or PRs
- Include references to issues if applicable (e.g., `Fixes #42`)
- You may use `[skip ci]` in a commit message to skip tests if you’re **only updating documentation**

---

## 🛡️ Branch Protection

The `master` branch is protected:

- ✅ Requires PRs for all changes
- ✅ Requires passing CI (`Build-and-test`)
- ✅ Requires branches to be up-to-date

You cannot push directly to `master` or bypass tests.

---

## 🧑‍💻 Useful Commands

Build from source:
```bash
mkdir build && cd build
cmake .. -DBUILD_EXAMPLES=ON -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

Run an example manually:
```bash
./build/examples/twoburn/twoburn
```

---

## 🙋 Need Help?

- Open a [GitHub Issue](https://github.com/PSOPT/psopt/issues)
- Contact the project maintainer

Thank you for helping improve PSOPT!
