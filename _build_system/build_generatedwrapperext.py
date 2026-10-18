"""Compiler-specific warning controls for generated SWIG extensions."""

from __future__ import annotations

from pathlib import Path
import subprocess

from setuptools.command.build_ext import build_ext


def deprecated_declaration_flag(compiler_type: str, command: list[str], banner: str) -> str | None:
    # Return the native option for deprecated-declaration warnings, if known.
    
    if compiler_type == "msvc":
        return "/wd4996"
    command_name = Path(command[0]).name.lower() if command else ""
    description = banner.lower()
    if "intel" in description or any(
        name in command_name for name in ("icc", "icpc")
    ):
        if "oneapi" not in description and "llvm" not in description and not any(
            name in command_name for name in ("icx", "icpx")
        ):
            # Classic Intel C++ warning: function was declared deprecated.
            return "-diag-disable=1478"
        return "-Wno-deprecated-declarations"
    if any(marker in description for marker in (
        "clang", "gcc", "g++", "gnu", "free software foundation",
    )):
        return "-Wno-deprecated-declarations"
    return None


class Build_NoDeprecationWarning(build_ext):
    # Apply known warning controls after setuptools selects a compiler.

    def build_extensions(self):
        command = list(getattr(self.compiler, "compiler_cxx", ()) or ())
        try:
            completed = subprocess.run(
                [*command, "--version"], capture_output=True, text=True,
                check=False,
            ) if command else None
        except OSError:
            completed = None
        banner = ((completed.stdout + completed.stderr)
                  if completed is not None else "")
        flag = deprecated_declaration_flag(
            self.compiler.compiler_type, command, banner
        )
        if flag is not None:
            for extension in self.extensions:
                if flag not in extension.extra_compile_args:
                    extension.extra_compile_args.append(flag)
        super().build_extensions()


def state_optimization_flag(compiler_type, banner):
    if compiler_type == 'unix' and any(marker in banner.lower() for marker in
            ('gcc', 'g++', 'clang', 'free software foundation')):
        return '-O3'
    return None


class Build_StateO3(Build_NoDeprecationWarning):
    def build_extensions(self):
        sources = frozenset(Path(source).resolve()
                            for extension in self.extensions
                            for source in getattr(extension, 'nsb_state_sources', ()))
        if not sources or self.compiler.compiler_type != 'unix':
            return super().build_extensions()

        command = list(getattr(self.compiler, 'compiler_cxx', ()) or ())
        try:
            completed = subprocess.run([*command, '--version'], capture_output=True,
                                       text=True, check=False) if command else None
        except OSError:
            completed = None
        banner = completed.stdout + completed.stderr if completed is not None else ''
        flag = state_optimization_flag(self.compiler.compiler_type, banner)
        if flag is None:
            return super().build_extensions()

        # UnixCompiler has no public per-source option API. Install one adapter
        # for the whole build, including parallel builds, and always restore it.
        original_compile = self.compiler._compile

        def compile_source(obj, src, ext, cc_args, extra_postargs, pp_opts):
            if Path(src).resolve() in sources:
                # The final optimization option overrides Python's usual -O2.
                # Copy options: setuptools can share them between sources.
                extra_postargs = [*(extra_postargs or ()), flag]
            return original_compile(obj, src, ext, cc_args, extra_postargs, pp_opts)

        self.compiler._compile = compile_source
        try:
            return super().build_extensions()
        finally:
            self.compiler._compile = original_compile
