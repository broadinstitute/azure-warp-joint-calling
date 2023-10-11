version 1.0

workflow GatherTest {
    input {
        Array[String] html_links
        String SAS_token
    }

    scatter(i in range(length(html_links))) {
        String inputs_concat = html_links[i] + SAS_token
    }

    scatter(i in range(50)){
        call GatherVcfs {
            input:
                input_vcfs = inputs_concat,
                disk_size = 50
        }
    }
}



task GatherVcfs {

  input {
    Array[String] input_vcfs
    String output_vcf_name = "test.vcf.gz"
    Int disk_size
    String gatk_docker = "mshand/genomesinthecloud:gatk_4.2.6.1"
  }

  command <<<
    set -euo pipefail

    # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
    # This argument disables expensive checks that the file headers contain the same set of
    # genotyped samples and that files are in order by position of first record.
    gatk --java-options "-Xms6000m -Xmx6500m" \
      GatherVcfsCloud \
      --ignore-safety-checks \
      --gather-type BLOCK \
      --input '~{sep="' --input '" input_vcfs}' \
      --output ~{output_vcf_name}

    tabix ~{output_vcf_name}
  >>>

  runtime {
    memory: "10000 MiB"
    cpu: "1"
    disk: disk_size + " GB"
    docker: gatk_docker
    maxRetries: 2
  }

  output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
  }
}