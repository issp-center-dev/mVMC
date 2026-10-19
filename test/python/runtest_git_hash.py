from __future__ import print_function

import os
import shutil
import subprocess
import sys
import tempfile

HASH = "0123456789abcdef0123456789abcdef01234567"


def run_script(cmake, script, source_dir, output):
    subprocess.run(
        [cmake, "-DSOURCE_DIR=" + source_dir, "-DOUTPUT=" + output,
         "-P", script],
        check=True)
    with open(output) as f:
        return f.read()


def git(source_dir, *arguments):
    proc = subprocess.run(
        ["git", "-c", "user.name=test", "-c", "user.email=test@example.com",
         "-c", "commit.gpgsign=false"] + list(arguments),
        cwd=source_dir, check=True, stdout=subprocess.PIPE, text=True)
    return proc.stdout.strip()


def write(filename, text):
    dirname = os.path.dirname(filename)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    with open(filename, "w") as f:
        f.write(text)


def main():
    if len(sys.argv) != 3:
        print("usage: {} <cmake> <git_hash.cmake>".format(sys.argv[0]))
        return -1
    cmake, script = sys.argv[1:3]

    status = 0
    workdir = tempfile.mkdtemp()

    def check(name, source_dir, expected):
        output = os.path.join(workdir, name + ".h")
        obtained = run_script(cmake, script, source_dir, output)
        if obtained == '#define MVMC_GIT_HASH "{}"\n'.format(expected):
            print("OK: {}: {!r}".format(name, expected))
            return 0
        print("ERROR: {}: expected {!r}, obtained {!r}".format(
            name, expected, obtained))
        return -1

    try:
        # copies which are neither a git repository nor an archive
        plain = os.path.join(workdir, "plain")
        os.makedirs(plain)
        status |= check("no_information", plain, "")

        unfilled = os.path.join(workdir, "unfilled")
        write(os.path.join(unfilled, "cmake", "git_archive.txt"),
              "$Format:%H$\n")
        status |= check("archive_not_filled_in", unfilled, "")

        # archive
        archive = os.path.join(workdir, "archive")
        write(os.path.join(archive, "cmake", "git_archive.txt"), HASH + "\n")
        status |= check("archive", archive, HASH[:8])

        # The header is not rewritten when the hash is the same.
        output = os.path.join(workdir, "archive.h")
        os.utime(output, (1000000000, 1000000000))
        run_script(cmake, script, archive, output)
        if os.stat(output).st_mtime == 1000000000:
            print("OK: header_not_rewritten")
        else:
            print("ERROR: the header was rewritten with the same hash")
            status = -1

        # git repository
        if shutil.which("git") is None:
            print("git is not found: the cases of a git repository are skipped")
            return status
        repository = os.path.join(workdir, "repository")
        write(os.path.join(repository, "tracked.txt"), "a\n")
        # an archive file in a repository is not read
        write(os.path.join(repository, "cmake", "git_archive.txt"),
              HASH + "\n")
        git(repository, "init", "-q")
        git(repository, "add", "tracked.txt", "cmake/git_archive.txt")
        git(repository, "commit", "-q", "-m", "test")
        head = git(repository, "rev-parse", "HEAD")
        status |= check("repository", repository, head[:8])

        write(os.path.join(repository, "untracked.txt"), "a\n")
        status |= check("repository_untracked", repository, head[:8])

        write(os.path.join(repository, "tracked.txt"), "b\n")
        status |= check("repository_dirty", repository, head[:8] + "-dirty")

        git(repository, "commit", "-q", "-a", "-m", "test 2")
        head = git(repository, "rev-parse", "HEAD")
        status |= check("repository_next_commit", repository, head[:8])
    finally:
        shutil.rmtree(workdir)
    return status


if __name__ == "__main__":
    sys.exit(main())
