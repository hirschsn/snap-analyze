// For license details see LICENSE.

#pragma once

#include <cstdio>
#include <string>
#include <cstdlib>
#include <cassert>

struct cfile {
    cfile() = delete;
    cfile(const cfile &) = delete;
    cfile &operator=(const cfile &) = delete;

    cfile(const std::string &fn, const char *mode)
        : f(fopen(fn.c_str(), mode)) {
        if (!f) {
            perror("open");
            exit(1);
        }
    }
    ~cfile() {
        assert(f); // No default constructor
        fclose(f);
    }
    operator FILE *() const { return f; }
    operator bool() const { return f != nullptr; }

  private:
    FILE *f = nullptr;
};
