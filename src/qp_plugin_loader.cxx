/*********************************************************************************************

This file is part of the PSOPT library, a software tool for computational optimal control

Copyright (C) 2009-2026 Victor M. Becerra

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA,
or visit http://www.gnu.org/licenses/

Author:    Professor Victor M. Becerra
Address:   University of Portsmouth
           School of Electrical and Mechanical Engineering
           Portsmouth PO1 3DJ
           United Kingdom
e-mail:    vmbecerra@vmb1.com

**********************************************************************************************/

//  Loading QP backend plugins.
//
//  The point of the exercise is the flag: RTLD_LOCAL. A plugin's symbols are not added
//  to the process's global namespace, so the ordering and factorisation code each one
//  carries stays its own, and two backends that would collide when linked together do
//  not collide when loaded side by side. Hidden visibility inside the plugin does the
//  other half of the job: its internal references are bound at link time to its own
//  copies rather than left for the dynamic linker to satisfy from whatever it finds
//  first. Neither half is sufficient alone.
//
//  Plugins are looked for in this order:
//
//    1. the directory named by PSOPT_QP_PLUGIN_PATH, if it is set, which is how a build
//       tree is tested before it is installed;
//    2. the directory the PSOPT library itself was loaded from, plus "psopt", which is
//       where they are installed and keeps a relocated installation working;
//    3. the install prefix compiled in at build time;
//    4. the build tree, so that tests and examples run before anything is installed;
//    5. the dynamic linker's own search path, by passing the bare file name.
//
//  A handle, once opened, is kept for the life of the process: a plugin has no state
//  between calls, and closing and reopening one for every subproblem would be absurd.

#include "psopt.h"
#include "psopt_qp_plugin.h"

#include <dlfcn.h>
#include <sys/stat.h>
#include <cctype>
#include <map>
#include <string>
#include <vector>

#ifndef PSOPT_QP_PLUGIN_INSTALL_DIR
#define PSOPT_QP_PLUGIN_INSTALL_DIR ""
#endif

#ifndef PSOPT_QP_PLUGIN_BUILD_DIR
#define PSOPT_QP_PLUGIN_BUILD_DIR ""
#endif

namespace {

struct LoadedPlugin {
    void*       handle = NULL;
    const char* (*name)(void) = NULL;
    int         (*abi)(void)  = NULL;
    int         (*solve)(const psopt_qp_problem*, psopt_qp_solution*) = NULL;
    int         (*environment_ok)(void) = NULL;   // optional
    std::string path;
    std::string error;
};

// The directory this library was loaded from, which is where the plugins sit beside it
// in an installed tree. dladdr resolves an address back to the object that provides it,
// so the library asks the dynamic linker where it is rather than being told at build
// time and then moved.
std::string library_directory()
{
    Dl_info info;
    if (dladdr((void*) &library_directory, &info) && info.dli_fname != NULL) {
        const std::string full(info.dli_fname);
        const size_t cut = full.find_last_of('/');
        if (cut != std::string::npos) return full.substr(0, cut);
    }
    return std::string();
}

std::vector<std::string> candidate_paths(const std::string& file)
{
    std::vector<std::string> out;

    const char* env = getenv("PSOPT_QP_PLUGIN_PATH");
    if (env != NULL && env[0] != '\0') out.push_back(std::string(env) + "/" + file);

    const std::string here = library_directory();
    if (!here.empty()) {
        out.push_back(here + "/" + file);
        out.push_back(here + "/psopt/" + file);
    }

    const std::string installed(PSOPT_QP_PLUGIN_INSTALL_DIR);
    if (!installed.empty()) out.push_back(installed + "/" + file);

    // The build tree, so that tests and examples run before anything is installed.
    const std::string built(PSOPT_QP_PLUGIN_BUILD_DIR);
    if (!built.empty()) out.push_back(built + "/" + file);

    out.push_back(file);          // leave it to the dynamic linker's search path
    return out;
}

LoadedPlugin& open_plugin(const std::string& backend)
{
    static std::map<std::string, LoadedPlugin> cache;

    std::map<std::string, LoadedPlugin>::iterator it = cache.find(backend);
    if (it != cache.end()) return it->second;

    LoadedPlugin p;

    // The option is spelled as the backend spells itself -- "ProxQP", "QPALM", "OSQP"
    // -- and the file is named in lower case, as libraries are.
    std::string tag = backend;
    for (size_t k = 0; k < tag.size(); k++)
        tag[k] = (char) tolower((unsigned char) tag[k]);
    const std::string file = "libpsopt_qp_" + tag + ".so";

    // Two failures wear the same face here and must not. A plugin that was never built
    // is absent; a plugin that was built but cannot be loaded -- because a shared library
    // it needs is not on the loader's path -- is present and broken. Only the second
    // names the real cause, in dlerror's text, and only the second is what a user with a
    // GALAHAD build and no LD_LIBRARY_PATH will meet.
    //
    // Reporting the LAST dlerror confuses them permanently. The last candidate tried is
    // the bare file name, left to the loader's own search, and its failure is always
    // "cannot open shared object file: No such file or directory" whatever went wrong
    // earlier. That message was once reported for a plugin sitting in the build tree
    // whose only problem was a missing libgalahad_double.so, and it sent the reader
    // looking for a file that was already there.
    //
    // So each candidate's error is kept as it happens, and the one from a candidate that
    // EXISTS is preferred over the rest.
    const std::vector<std::string> paths = candidate_paths(file);
    std::string existed_but_failed;
    std::string looked_in;

    for (size_t k = 0; k < paths.size() && p.handle == NULL; k++) {
        // RTLD_LOCAL is the whole point; RTLD_NOW so that a plugin missing a symbol
        // says so on load rather than in the middle of a solve.
        p.handle = dlopen(paths[k].c_str(), RTLD_NOW | RTLD_LOCAL);
        if (p.handle != NULL) { p.path = paths[k]; break; }

        // dlerror() must be read now: it is cleared by the next call, so an error kept
        // until after the loop is not the error of the attempt one means.
        const char* why = dlerror();

        if (!looked_in.empty()) looked_in += ", ";
        looked_in += paths[k];

        struct stat st;
        if (existed_but_failed.empty() && stat(paths[k].c_str(), &st) == 0 && why != NULL)
            existed_but_failed = paths[k] + ": " + why;
    }

    if (p.handle == NULL) {
        if (!existed_but_failed.empty()) {
            p.error = "could not load the QP backend plugin " + file
                    + ", which is present but will not open (" + existed_but_failed
                    + "). The parenthesis is the loader's own words and is the thing to "
                      "read: a missing shared library means the backend is not on the "
                      "loader's path, and an executable-stack refusal means the backend "
                      "was built asking for one, which glibc no longer grants at dlopen "
                      "and which the backend has to be rebuilt without.";
        }
        else {
            p.error = "could not load the QP backend plugin " + file
                    + ", which was not found. Looked in: " + looked_in
                    + ". Set PSOPT_QP_PLUGIN_PATH to the directory containing it, or "
                      "build PSOPT with that backend enabled.";
        }
    }
    else {
        p.abi   = (int (*)(void)) dlsym(p.handle, "psopt_qp_abi_version");
        p.name  = (const char* (*)(void)) dlsym(p.handle, "psopt_qp_name");
        p.solve = (int (*)(const psopt_qp_problem*, psopt_qp_solution*))
                      dlsym(p.handle, "psopt_qp_solve");

        // Optional, and provided only by a backend with a requirement it cannot satisfy
        // for itself. GALAHAD's default linear solver needs OMP_CANCELLATION and
        // OMP_PROC_BIND set in the environment, which the OpenMP runtime reads when it
        // first initialises -- long before a plugin is loaded, so nothing the plugin
        // does can put it right. Unmet, it does not announce itself: the solver returns
        // an iteration-limit code and a vector that is not a solution, on problems of
        // two variables. Asking here turns a day's diagnosis into a sentence.
        p.environment_ok = (int (*)(void)) dlsym(p.handle, "psopt_qp_environment_ok");

        if (p.abi == NULL || p.name == NULL || p.solve == NULL)
            p.error = "the QP backend plugin " + p.path
                    + " does not provide the expected entry points";
        else if (p.abi() != PSOPT_QP_ABI_VERSION)
            p.error = "the QP backend plugin " + p.path + " was built against plugin ABI "
                    + std::to_string(p.abi()) + ", and this PSOPT expects "
                    + std::to_string(PSOPT_QP_ABI_VERSION);
        else if (p.environment_ok != NULL && p.environment_ok() == 0)
            p.error = std::string("the ") + p.name() + " backend cannot run in this "
                      "environment. Its linear solver requires OMP_CANCELLATION=TRUE and "
                      "OMP_PROC_BIND=TRUE, and they must be set before the program starts, "
                      "because the OpenMP runtime reads them once when it initialises. "
                      "Set both in the environment and run again.";
    }

    cache[backend] = p;
    // The map holds the entry by value; return the stored one, not the local.
    return cache[backend];
}

} // anonymous namespace


// Solve a subproblem through a named backend. Returns false and leaves a message in
// 'message' if the plugin could not be loaded; a subproblem the backend declines is
// reported through solution->status and is not an error here.
bool psopt_qp_plugin_solve(const std::string& backend,
                           const psopt_qp_problem* problem,
                           psopt_qp_solution*      solution,
                           std::string&            message)
{
    LoadedPlugin& p = open_plugin(backend);

    if (!p.error.empty()) { message = p.error; return false; }

    p.solve(problem, solution);
    return true;
}

// Whether a backend is available, without solving anything. Used by validate() so that
// a missing plugin is reported when the problem is set up rather than at the first
// subproblem, and by the tests.
bool psopt_qp_plugin_available(const std::string& backend, std::string& message)
{
    LoadedPlugin& p = open_plugin(backend);
    if (!p.error.empty()) { message = p.error; return false; }
    message = p.name();
    return true;
}
