version 1.0

workflow CreateSampleMapFile{
    
    input {
        Array[String] sample_names
        Array[String] gvcfs
    }

    call CreateSampleMap {
        input:
            sample_names=sample_names,
            input_gvcfs=gvcfs
    }
    output {
        File sample_map_file = CreateSampleMap.output_map
    }
}


task CreateSampleMap {
    input {
        Array[String] sample_names
        Array[String] input_gvcfs
    }

    command <<<

        set -xe

        python << CODE
        import re

        gvcfs = ['~{sep="','" input_gvcfs}']
        samples = ['~{sep="','" sample_names}']

        if len(gvcfs) != len(samples):
          print(f"error! len(gvcfs)={len(gvcfs)} which is different from len(sample_names)={len(samples)}. Quiting!")    
          exit(1)

        with open("inputs.list", "w") as fi:
          for i in range(len(gvcfs)):
            container = re.search("sc-.{8}-.{4}-.{4}-.{4}-.{12}", gvcfs[i])
            account = re.search("(lz.*).blob", gvcfs[i])
            gvcf_path = re.search("(sc-.{8}-.{4}-.{4}-.{4}-.{12})(\/.*)", gvcfs[i])
            az_path = "az://" + container.group(0) + "@" + account.group(1) + ".blob.core.windows.net" + gvcf_path.group(2)
            fi.write(samples[i] + "\t" + az_path + "\n") 

        CODE

    >>>
    runtime {
        docker: "python:3.6"
        memory: "7 GB"
        cpu: "2"
        disks: "local-disk " + 50 + " HDD"
    }
    output {
        File output_map = "inputs.list"
    }
}