# click_mergevcfs

[![pypi badge][pypi_badge]][pypi_base]
[![travis badge][travis_badge]][travis_base]
[![codecov badge][codecov_badge]][codecov_base]
[![docker badge][docker_badge]][docker_base]
[![docker badge][automated_badge]][docker_base]

Merge vcfs files from multiple different callers.

## Installation

        pip install click_mergevcfs

## Usage

        singularity run [SINGULARITY-OPTIONS] <image> [PIPELINE-OPTIONS]
        or
        singularity run [SINGULARITY-OPTIONS] docker://leukgen/click_mergevcfs [PIPELINE-OPTIONS]

## Example

        singularity \
                run \
                --containall \
                --bind /ifs \
                {singularity image} \
                --vcf {path to vcf1} \
                --vcf {path to vcf2} \
                --vcf {path to vcf3} \
                --out {path to merged output vcf.gz} \
                --snv

## Options

| Option      | Description                                   |
| ----------- | --------------------------------------------- |
| --vcf       | Path to a input vcf file                      |
| --out       | Path to the output file                       |
| --snv       | If the input vcf files contain only snvs      |
| --indel     | If the input vcf files contain only indels    |
| --sv        | If the input vcf files contain only svs       |
| --reference | Path to the reference fasta                   |
| --temp      | Temporary working directory                   |

The following options are required to run caveman postprocessing:

| Option            | Description                                   |
| ----------- | --------------------------------------------- |
| --caveman_flagged_out | Path to caveman flagged output file |
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