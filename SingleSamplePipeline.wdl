version 1.0

workflow SingleSamplePipeline {
    input {
        String name
        Int value
    }

    call task1 {
        input: 
            name = name, 
            value = value
    }

    output {
        File out = task1.out
    }
}

task task1 {
    input {
        String name
        Int value
    }

    command {
        echo "Hello, ${name}! Your value is ${value}."
    }

    output {
        File out = stdout()
    }

    runtime {
        docker: "ubuntu:latest"
    }
}