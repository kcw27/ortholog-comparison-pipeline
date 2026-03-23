process foo {
    output:
    val "foo"

    script:
    """
    echo "foo"
    """
}

process bar {
    output:
    val "bar"

    script:
    """
    echo "bar"
    """
}

process handleGenomeFlag {
    input:
    val flag
    
    script:
    """
    if [ "$flag" == "true" ]; then
        echo "Processing with genome IDs"
        # Your commands here
    else
        echo "Processing without genome IDs"
        # Your commands here
    fi
    """
}