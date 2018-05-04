# click_mergevcfs

[![pypi badge][pypi_badge]][pypi_base]
[![travis badge][travis_badge]][travis_base]
[![codecov badge][codecov_badge]][codecov_base]
[![docker badge][docker_badge]][docker_base]
[![docker badge][automated_badge]][docker_base]

Merge vcfs files from multiple different callers.

## Installation

        pip install click_mergevcfs

## Example

This tool is designed to be used with containers. Run click_mergevcfs with a container:

        docker run --volumes local/path:/data leukgen/click_mergevcfs \
                --vcf /path/to/caller1.snvs.vcf \
                --vcf /path/to/caller2.snvs.vcf \
                --vcf /path/to/caller3.snvs.vcf \
                --out /path/to/merged_output.vcf \
                --reference /path/to/reference.fasta \
                --snv \

## Options

| Option      | Description                                   |
| ----------- | --------------------------------------------- |
| --vcf       | Path to a input vcf file                      |
| --outdir    | Path to the output file                       |
| --snv       | If the input vcf files contain snvs or indels |
| --sv        | If the input vcf files contain svs            |
| --reference | Path to the reference fasta                   |
| --no_flag   | Don't apply custom postprocessing             |
| --temp      | Temporary working directory                   |

The following options are required to run caveman postprocessing:

| Option            | Description                                   |
| ----------- | --------------------------------------------- |
| --normal_bam      | Path to normal bam |
| --tumor_bam       | Path to tumor bam |
| --bedFileLoc      | Path to a folder containing centromeric, snp, hi seq depth, simple repeat bed files |
| --indelBed        | A bed file containing germline indels to filter on |
| --unmatchedVCFLoc | Path to folder containing unmatched VCF PON |
| --annoBedLoc      | Path to bed files containing annotatable regions and coding regions |

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
[automated_badge]: https://img.shields.io/docker/automated/leukgen/click_mergevcfs.svg
[docker_base]: https://hub.docker.com/r/leukgen/click_mergevcfs
[docker_badge]: https://img.shields.io/docker/build/leukgen/click_mergevcfs.svg
[pypi_badge]: https://img.shields.io/pypi/v/click_mergevcfs.svg
[pypi_base]: https://pypi.python.org/pypi/click_mergevcfs
[travis_badge]: https://img.shields.io/travis/leukgen/click_mergevcfs.svg
[travis_base]: https://travis-ci.org/leukgen/click_mergevcfs