# click_mergevcfs

[![pypi badge][pypi_badge]][pypi_base]
[![travis badge][travis_badge]][travis_base]
[![codecov badge][codecov_badge]][codecov_base]
[![docker badge][docker_badge]][docker_base]
[![docker badge][automated_badge]][docker_base]
[![code formatting][black_badge]][black_base]

Merge vcfs files from multiple different callers.

Due to different caller output vcf files have with slightly different formats, merging vcfs is not trivial; several steps need to be done before merging vcfs:

- First, the order of the last 2 columns in the output vcf of each callers must be NOMRAL TUMOR; if that's not the case (ex. mutect), click_mergevcfs will rearrange the column order.
- Second, multiallelic variants are formatted differently by each caller. For example, strelka's format is one REF, multiple ALT; mutect's format is one REF, one ALT, multiple records. click_mergevcfs will break mutliple ALT into multiple records each with one ALT.
- Then, click_mergevcfs will add PASSED_{caller} information under the INFO field if a record has the flag PASS. This is to distinguish passed variants once variants are merged from callers.
- Finally, vcf files are merged using `vcf-merge --collapse none vcf1 vcf2 vcf3`. This will ensure the output vcf keeps all conflicting calls. Headers are also formatted to include information from all callers.
- For SNVs only, we apply caveman flagging to all variants from all callers. This is done by calling an altered `cgpFlagCaVEMan.pl`, which was changed to accommmedate variants that are called by mutect/strelka that have zero supporting reads in the original bam because mutect/strelka does on-the-fly indel realginment. When flagging WGS vcfs, caveman flagging is implementated in a parallelized and memory-efficient manner to achieve reasonable runtime. It is done by spliting the input merged vcf into many smaller vcf, apply `cgpFlagCaVEMan.pl` in parallel, and concat output vcfs into one output vcf. Additional safety measures were implemented so that it only success if the number of variants in the output vcf is the same as that in the input vcf. All intermediate files are deleted.

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
                --scratch {recommend using /tmp} \
                --workdir {recommend using $TMP_DIR} \
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
[black_badge]: https://img.shields.io/badge/code%20style-black-000000.svg
[black_base]: https://github.com/ambv/black
