# click_mergevcfs

[![pypi badge][pypi_badge]][pypi_base]
[![travis badge][travis_badge]][travis_base]
[![codecov badge][codecov_badge]][codecov_base]

Merge multiple vcf files and apply custom post-processing

## Features

* 📦 &nbsp; **Easy Installation**

        pip install click_mergevcfs

* 🍉 &nbsp; **Usage Documentation**

        click_mergevcfs --help

* 🐳 &nbsp; **Containers Support**

        # docker usage
        docker run --volume /shared_fs:/shared_fs --interactive --tty \
            click_mergevcfs-image
            [click_mergevcfs options]

        # singularity usage
        singularity run --workdir /shared_fs/tmp --bind /shared_fs:/shared_fs \
            click_mergevcfs-singularity-image-path
            [click_mergevcfs options]

## Contributing

Contributions are welcome, and they are greatly appreciated, check our [contributing guidelines](.github/CONTRIBUTING.md)!

## Credits

This package was created using [Cookiecutter] and the
[leukgen/cookiecutter-toil] project template.

<!-- References -->
[singularity]: http://singularity.lbl.gov/
[docker2singularity]: https://github.com/singularityware/docker2singularity
[cookiecutter]: https://github.com/audreyr/cookiecutter
[leukgen/cookiecutter-toil]: https://github.com/leukgen/cookiecutter-toil

<!-- Badges -->
[codecov_badge]: https://codecov.io/gh/leukgen/click_mergevcfs/branch/master/graph/badge.svg
[codecov_base]: https://codecov.io/gh/leukgen/click_mergevcfs
[pypi_badge]: https://img.shields.io/pypi/v/click_mergevcfs.svg
[pypi_base]: https://pypi.python.org/pypi/click_mergevcfs
[travis_badge]: https://img.shields.io/travis/leukgen/click_mergevcfs.svg
[travis_base]: https://travis-ci.org/leukgen/click_mergevcfs
