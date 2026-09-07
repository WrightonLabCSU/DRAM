import java.nio.file.Files
import java.nio.file.Path

/*
 * Return the total number of bytes represented by a path, directory, or nested
 * collection of paths. Directory contents are summed recursively.
 */
def resourceBytes(value) {
    if (value == null) {
        return 0L
    }
    if (value instanceof Map) {
        return value.values().collect { resourceBytes(it) }.sum(0L) as long
    }
    if (value instanceof Collection) {
        return value.collect { resourceBytes(it) }.sum(0L) as long
    }
    if (value.getClass().isArray()) {
        return value.collect { resourceBytes(it) }.sum(0L) as long
    }
    if (!(value instanceof Path)) {
        throw new IllegalArgumentException("Cannot calculate resource size for ${value.getClass().name}")
    }
    if (!Files.isDirectory(value)) {
        return Files.size(value)
    }

    def paths = Files.walk(value)
    try {
        return paths
            .filter { Files.isRegularFile(it) }
            .mapToLong { Files.size(it) }
            .sum()
    }
    finally {
        paths.close()
    }
}

/*
 * Resource classes intentionally use binary units because Nextflow's GiB-sized
 * resource values are also binary quantities.
 */
def resourceClass(long bytes) {
    final long ONE_GIB = 1024L * 1024L * 1024L
    final long TWENTY_GIB = 20L * ONE_GIB

    if (bytes <= ONE_GIB) {
        return 'small'
    }
    if (bytes <= TWENTY_GIB) {
        return 'medium'
    }
    return 'large'
}
