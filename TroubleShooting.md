# Why is my SAM file empty after alignment? 
## Checklist 
- ✅ Confirm the FASTQ files actually contain reads 
    ```
    $ zcat SRR062634_1.filt.fastq.gz | wc -l
    96595972
    ```
    ```
    $ zcat SRR062634_2.filt.fastq.gz | wc -l
    96595972
    ```
- ✅ Check the files exist & paths are correct
- Make sure the reference index exists and ${ref} points to the indexed reference basename
    ```
    $ ls -lh /home/dmatute/Projects/Variant_Calling/data/hg38.fa.* 2>/dev/null || ls -lh /home/dmatute/Projects/Variant_Calling/data/hg38.fa.bwt /home/dmatute/Projects/Variant_Calli
    ng/data/hg38.fa.pac /home/dmatute/Projects/Variant_Calling/data/hg38.fa.ann ${ref}.amb /home/dmatute/Projects/Variant_Calling/data/hg38.fa.sa 2>/dev/null
    -rw-r--r-- 1 dmatute dmatute  21K Oct 10 18:53 /home/dmatute/Projects/Variant_Calling/data/hg38.fa.amb
    -rw-r--r-- 1 dmatute dmatute  22K Oct 10 18:53 /home/dmatute/Projects/Variant_Calling/data/hg38.fa.ann
    -rw-r--r-- 1 dmatute dmatute 3.0G Oct 10 18:53 /home/dmatute/Projects/Variant_Calling/data/hg38.fa.bwt
    -rw-r--r-- 1 dmatute dmatute  19K Oct  6 21:57 /home/dmatute/Projects/Variant_Calling/data/hg38.fa.fai
    -rw-r--r-- 1 dmatute dmatute 766M Oct 10 18:53 /home/dmatute/Projects/Variant_Calling/data/hg38.fa.pac
    -rw-r--r-- 1 dmatute dmatute 1.5G Oct 10 19:30 /home/dmatute/Projects/Variant_Calling/data/hg38.fa.sa
    (base) dmatute@DanyDell:~/Projects/Variant_Calling$ 
    ```
- ✅ Run a tiny test alignment (stdout + stderr) — do not redirect to file for this test
    <details>

        (gatk-env) dmatute@DanyDell:~/Projects/Variant_Calling/test$ bwa mem -t 2 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" /home/dmatute/Projects/Variant_Calling/data/hg38.fa ../reads/SRR062634_1.filt.fastq.gz ../reads/SRR062634_2.filt.fastq.gz | head -n 40
        [M::bwa_idx_load_from_disk] read 0 ALT contigs
        @HD     VN:1.5  SO:unsorted     GO:query
        @SQ     SN:chr1 LN:248956422
        @SQ     SN:chr10        LN:133797422
        @SQ     SN:chr11        LN:135086622
        @SQ     SN:chr11_KI270721v1_random      LN:100316
        @SQ     SN:chr12        LN:133275309
        @SQ     SN:chr13        LN:114364328
        @SQ     SN:chr14        LN:107043718
        @SQ     SN:chr14_GL000009v2_random      LN:201709
        @SQ     SN:chr14_GL000225v1_random      LN:211173
        @SQ     SN:chr14_KI270722v1_random      LN:194050
        @SQ     SN:chr14_GL000194v1_random      LN:191469
        @SQ     SN:chr14_KI270723v1_random      LN:38115
        @SQ     SN:chr14_KI270724v1_random      LN:39555
        @SQ     SN:chr14_KI270725v1_random      LN:172810
        @SQ     SN:chr14_KI270726v1_random      LN:43739
        @SQ     SN:chr15        LN:101991189
        @SQ     SN:chr15_KI270727v1_random      LN:448248
        @SQ     SN:chr16        LN:90338345
        @SQ     SN:chr16_KI270728v1_random      LN:1872759
        @SQ     SN:chr17        LN:83257441
        @SQ     SN:chr17_GL000205v2_random      LN:185591
        @SQ     SN:chr17_KI270729v1_random      LN:280839
        @SQ     SN:chr17_KI270730v1_random      LN:112551
        @SQ     SN:chr18        LN:80373285
        @SQ     SN:chr19        LN:58617616
        @SQ     SN:chr1_KI270706v1_random       LN:175055
        @SQ     SN:chr1_KI270707v1_random       LN:32032
        @SQ     SN:chr1_KI270708v1_random       LN:127682
        @SQ     SN:chr1_KI270709v1_random       LN:66860
        @SQ     SN:chr1_KI270710v1_random       LN:40176
        @SQ     SN:chr1_KI270711v1_random       LN:42210
        @SQ     SN:chr1_KI270712v1_random       LN:176043
        @SQ     SN:chr1_KI270713v1_random       LN:40745
        @SQ     SN:chr1_KI270714v1_random       LN:41717
        @SQ     SN:chr2 LN:242193529
        @SQ     SN:chr20        LN:64444167
        @SQ     SN:chr21        LN:46709983
        @SQ     SN:chr22        LN:50818468
        @SQ     SN:chr22_KI270731v1_random      LN:150754
        [M::process] read 200000 sequences (20000000 bp)...
        [M::process] read 200000 sequences (20000000 bp)...
        [M::mem_pestat] # candidate unique pairs for (FF, FR, RF, RR): (2, 80678, 15, 1)
        [M::mem_pestat] skip orientation FF as there are not enough pairs
        [M::mem_pestat] analyzing insert size distribution for orientation FR...
        [M::mem_pestat] (25, 50, 75) percentile: (168, 181, 194)
        [M::mem_pestat] low and high boundaries for computing mean and std.dev: (116, 246)
        [M::mem_pestat] mean and std.dev: (180.52, 19.47)
        [M::mem_pestat] low and high boundaries for proper pairs: (90, 272)
        [M::mem_pestat] analyzing insert size distribution for orientation RF...
        [M::mem_pestat] (25, 50, 75) percentile: (186, 202, 211)
        [M::mem_pestat] low and high boundaries for computing mean and std.dev: (136, 261)
        [M::mem_pestat] mean and std.dev: (197.69, 14.46)
        [M::mem_pestat] low and high boundaries for proper pairs: (111, 286)
        [M::mem_pestat] skip orientation RR as there are not enough pairs
        [M::mem_pestat] skip orientation RF
        [M::mem_process_seqs] Processed 200000 reads in 82.414 CPU sec, 43.762 real sec
        
    </details>

Everything seems fine. 👍
You just stopped the BWA alignment early. 

➡️ Re run BWA
```bwa mem -t 8 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" \
/home/dmatute/Projects/Variant_Calling/data/hg38.fa \
../reads/SRR062634_1.filt.fastq.gz \
../reads/SRR062634_2.filt.fastq.gz \
> ../aligned_reads/SRR062634.paired.sam
```

# Why is ```gatk MarkDuplicatesSpar``` not working

## This is the error
<details>
    ```

    $ gatk MarkDuplicatesSpark -I SRR062634_sorted.bam -O SRR062634_sorted_deduo_reads.baw

    Using GATK jar /home/dmatute/miniconda3/envs/gatk-env/share/gatk4-4.3.0.0-0/gatk-package-4.3.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /home/dmatute/miniconda3/envs/gatk-env/share/gatk4-4.3.0.0-0/gatk-package-4.3.0.0-local.jar MarkDuplicatesSpark -I SRR062634_sorted.bam -O SRR062634_sorted_deduo_reads.baw
    22:41:07.441 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/home/dmatute/miniconda3/envs/gatk-env/share/gatk4-4.3.0.0-0/gatk-package-4.3.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    22:41:27.773 INFO  MarkDuplicatesSpark - ------------------------------------------------------------
    22:41:27.774 INFO  MarkDuplicatesSpark - The Genome Analysis Toolkit (GATK) v4.3.0.0
    22:41:27.774 INFO  MarkDuplicatesSpark - For support and documentation go to https://software.broadinstitute.org/gatk/
    22:41:27.774 INFO  MarkDuplicatesSpark - Executing as dmatute@DanyDell on Linux v6.6.87.2-microsoft-standard-WSL2 amd64
    22:41:27.774 INFO  MarkDuplicatesSpark - Java runtime: OpenJDK 64-Bit Server VM v11.0.27+6
    22:41:27.775 INFO  MarkDuplicatesSpark - Start Date/Time: October 15, 2025 at 10:41:07 PM EDT
    22:41:27.775 INFO  MarkDuplicatesSpark - ------------------------------------------------------------
    22:41:27.775 INFO  MarkDuplicatesSpark - ------------------------------------------------------------
    22:41:27.775 INFO  MarkDuplicatesSpark - HTSJDK Version: 3.0.1
    22:41:27.776 INFO  MarkDuplicatesSpark - Picard Version: 2.27.5
    22:41:27.776 INFO  MarkDuplicatesSpark - Built for Spark Version: 2.4.5
    22:41:27.776 INFO  MarkDuplicatesSpark - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    22:41:27.776 INFO  MarkDuplicatesSpark - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    22:41:27.776 INFO  MarkDuplicatesSpark - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    22:41:27.776 INFO  MarkDuplicatesSpark - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    22:41:27.776 INFO  MarkDuplicatesSpark - Deflater: IntelDeflater
    22:41:27.776 INFO  MarkDuplicatesSpark - Inflater: IntelInflater
    22:41:27.776 INFO  MarkDuplicatesSpark - GCS max retries/reopens: 20
    22:41:27.777 INFO  MarkDuplicatesSpark - Requester pays: disabled
    22:41:27.777 INFO  MarkDuplicatesSpark - Initializing engine
    22:41:27.777 INFO  MarkDuplicatesSpark - Done initializing engine
    Using Spark's default log4j profile: org/apache/spark/log4j-defaults.properties
    25/10/15 22:41:28 WARN Utils: Your hostname, DanyDell resolves to a loopback address: 127.0.1.1; using 10.255.255.254 instead (on interface lo)
    25/10/15 22:41:28 WARN Utils: Set SPARK_LOCAL_IP if you need to bind to another address
    WARNING: An illegal reflective access operation has occurred
    WARNING: Illegal reflective access by org.apache.spark.unsafe.Platform (file:/home/dmatute/miniconda3/envs/gatk-env/share/gatk4-4.3.0.0-0/gatk-package-4.3.0.0-local.jar) to method java.nio.Bits.unaligned()
    WARNING: Please consider reporting this to the maintainers of org.apache.spark.unsafe.Platform
    WARNING: Use --illegal-access=warn to enable warnings of further illegal reflective access operations
    WARNING: All illegal access operations will be denied in a future release
    25/10/15 22:41:38 INFO SparkContext: Running Spark version 2.4.5
    25/10/15 22:41:38 WARN NativeCodeLoader: Unable to load native-hadoop library for your platform... using builtin-java classes where applicable
    25/10/15 22:41:38 INFO SparkContext: Submitted application: MarkDuplicatesSpark
    25/10/15 22:41:39 INFO SecurityManager: Changing view acls to: dmatute
    25/10/15 22:41:39 INFO SecurityManager: Changing modify acls to: dmatute
    25/10/15 22:41:39 INFO SecurityManager: Changing view acls groups to: 
    25/10/15 22:41:39 INFO SecurityManager: Changing modify acls groups to: 
    25/10/15 22:41:39 INFO SecurityManager: SecurityManager: authentication disabled; ui acls disabled; users  with view permissions: Set(dmatute); groups with view permissions: Set(); users  with modify permissions: Set(dmatute); groups with modify permissions: Set()
    25/10/15 22:41:39 INFO Utils: Successfully started service 'sparkDriver' on port 43845.
    25/10/15 22:41:39 INFO SparkEnv: Registering MapOutputTracker
    25/10/15 22:41:39 INFO SparkEnv: Registering BlockManagerMaster
    25/10/15 22:41:39 INFO BlockManagerMasterEndpoint: Using org.apache.spark.storage.DefaultTopologyMapper for getting topology information
    25/10/15 22:41:39 INFO BlockManagerMasterEndpoint: BlockManagerMasterEndpoint up
    25/10/15 22:41:39 INFO DiskBlockManager: Created local directory at /tmp/blockmgr-474834b9-fa98-433c-b7d2-9fa1b8ac3568
    25/10/15 22:41:39 INFO MemoryStore: MemoryStore started with capacity 988.8 MB
    25/10/15 22:41:39 INFO SparkEnv: Registering OutputCommitCoordinator
    25/10/15 22:41:40 INFO Utils: Successfully started service 'SparkUI' on port 4040.
    25/10/15 22:41:40 INFO SparkUI: Bound SparkUI to 0.0.0.0, and started at http://10.255.255.254:4040
    25/10/15 22:41:40 INFO Executor: Starting executor ID driver on host localhost
    25/10/15 22:41:40 INFO Utils: Successfully started service 'org.apache.spark.network.netty.NettyBlockTransferService' on port 40579.
    25/10/15 22:41:40 INFO NettyBlockTransferService: Server created on 10.255.255.254:40579
    25/10/15 22:41:40 INFO BlockManager: Using org.apache.spark.storage.RandomBlockReplicationPolicy for block replication policy
    25/10/15 22:41:40 INFO BlockManagerMaster: Registering BlockManager BlockManagerId(driver, 10.255.255.254, 40579, None)
    25/10/15 22:41:40 INFO BlockManagerMasterEndpoint: Registering block manager 10.255.255.254:40579 with 988.8 MB RAM, BlockManagerId(driver, 10.255.255.254, 40579, None)
    25/10/15 22:41:40 INFO BlockManagerMaster: Registered BlockManager BlockManagerId(driver, 10.255.255.254, 40579, None)
    25/10/15 22:41:40 INFO BlockManager: Initialized BlockManager: BlockManagerId(driver, 10.255.255.254, 40579, None)
    22:41:40.792 INFO  MarkDuplicatesSpark - Spark verbosity set to INFO (see --spark-verbosity argument)
    25/10/15 22:41:40 INFO GoogleHadoopFileSystemBase: GHFS version: 1.9.4-hadoop3
    25/10/15 22:41:41 INFO MemoryStore: Block broadcast_0 stored as values in memory (estimated size 172.9 KB, free 988.6 MB)
    25/10/15 22:41:42 INFO MemoryStore: Block broadcast_0_piece0 stored as bytes in memory (estimated size 35.5 KB, free 988.6 MB)
    25/10/15 22:41:42 INFO BlockManagerInfo: Added broadcast_0_piece0 in memory on 10.255.255.254:40579 (size: 35.5 KB, free: 988.8 MB)
    25/10/15 22:41:42 INFO SparkContext: Created broadcast 0 from newAPIHadoopFile at PathSplitSource.java:96
    25/10/15 22:41:42 INFO BlockManagerInfo: Removed broadcast_0_piece0 on 10.255.255.254:40579 in memory (size: 35.5 KB, free: 988.8 MB)
    25/10/15 22:41:42 INFO MemoryStore: Block broadcast_1 stored as values in memory (estimated size 172.9 KB, free 988.6 MB)
    25/10/15 22:41:42 INFO MemoryStore: Block broadcast_1_piece0 stored as bytes in memory (estimated size 35.5 KB, free 988.6 MB)
    25/10/15 22:41:42 INFO BlockManagerInfo: Added broadcast_1_piece0 in memory on 10.255.255.254:40579 (size: 35.5 KB, free: 988.8 MB)
    25/10/15 22:41:42 INFO SparkContext: Created broadcast 1 from newAPIHadoopFile at PathSplitSource.java:96
    25/10/15 22:41:42 INFO FileInputFormat: Total input files to process : 1
    25/10/15 22:41:42 INFO SparkUI: Stopped Spark web UI at http://10.255.255.254:4040
    25/10/15 22:41:42 INFO MapOutputTrackerMasterEndpoint: MapOutputTrackerMasterEndpoint stopped!
    25/10/15 22:41:42 INFO MemoryStore: MemoryStore cleared
    25/10/15 22:41:42 INFO BlockManager: BlockManager stopped
    25/10/15 22:41:42 INFO BlockManagerMaster: BlockManagerMaster stopped
    25/10/15 22:41:42 INFO OutputCommitCoordinator$OutputCommitCoordinatorEndpoint: OutputCommitCoordinator stopped!
    25/10/15 22:41:42 INFO SparkContext: Successfully stopped SparkContext
    22:41:42.752 INFO  MarkDuplicatesSpark - Shutting down engine
    [October 15, 2025 at 10:41:42 PM EDT] org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSpark done. Elapsed time: 0.59 minutes.
    Runtime.totalMemory()=154140672
    java.lang.IllegalArgumentException: Unsupported class file major version 55
            at org.apache.xbean.asm6.ClassReader.<init>(ClassReader.java:166)
            at org.apache.xbean.asm6.ClassReader.<init>(ClassReader.java:148)
            at org.apache.xbean.asm6.ClassReader.<init>(ClassReader.java:136)
            at org.apache.xbean.asm6.ClassReader.<init>(ClassReader.java:237)
            at org.apache.spark.util.ClosureCleaner$.getClassReader(ClosureCleaner.scala:49)
            at org.apache.spark.util.FieldAccessFinder$$anon$3$$anonfun$visitMethodInsn$2.apply(ClosureCleaner.scala:517)
            at org.apache.spark.util.FieldAccessFinder$$anon$3$$anonfun$visitMethodInsn$2.apply(ClosureCleaner.scala:500)
            at scala.collection.TraversableLike$WithFilter$$anonfun$foreach$1.apply(TraversableLike.scala:733)
            at scala.collection.mutable.HashMap$$anon$1$$anonfun$foreach$2.apply(HashMap.scala:134)
            at scala.collection.mutable.HashMap$$anon$1$$anonfun$foreach$2.apply(HashMap.scala:134)
            at scala.collection.mutable.HashTable$class.foreachEntry(HashTable.scala:236)
            at scala.collection.mutable.HashMap.foreachEntry(HashMap.scala:40)
            at scala.collection.mutable.HashMap$$anon$1.foreach(HashMap.scala:134)
            at scala.collection.TraversableLike$WithFilter.foreach(TraversableLike.scala:732)
            at org.apache.spark.util.FieldAccessFinder$$anon$3.visitMethodInsn(ClosureCleaner.scala:500)
            at org.apache.xbean.asm6.ClassReader.readCode(ClassReader.java:2175)
            at org.apache.xbean.asm6.ClassReader.readMethod(ClassReader.java:1238)
            at org.apache.xbean.asm6.ClassReader.accept(ClassReader.java:631)
            at org.apache.xbean.asm6.ClassReader.accept(ClassReader.java:355)
            at org.apache.spark.util.ClosureCleaner$$anonfun$org$apache$spark$util$ClosureCleaner$$clean$14.apply(ClosureCleaner.scala:307)
            at org.apache.spark.util.ClosureCleaner$$anonfun$org$apache$spark$util$ClosureCleaner$$clean$14.apply(ClosureCleaner.scala:306)
            at scala.collection.immutable.List.foreach(List.scala:392)
            at org.apache.spark.util.ClosureCleaner$.org$apache$spark$util$ClosureCleaner$$clean(ClosureCleaner.scala:306)
            at org.apache.spark.util.ClosureCleaner$.clean(ClosureCleaner.scala:162)
            at org.apache.spark.SparkContext.clean(SparkContext.scala:2326)
            at org.apache.spark.SparkContext.runJob(SparkContext.scala:2100)
            at org.apache.spark.SparkContext.runJob(SparkContext.scala:2126)
            at org.apache.spark.rdd.RDD$$anonfun$collect$1.apply(RDD.scala:990)
            at org.apache.spark.rdd.RDDOperationScope$.withScope(RDDOperationScope.scala:151)
            at org.apache.spark.rdd.RDDOperationScope$.withScope(RDDOperationScope.scala:112)
            at org.apache.spark.rdd.RDD.withScope(RDD.scala:385)
            at org.apache.spark.rdd.RDD.collect(RDD.scala:989)
            at org.apache.spark.RangePartitioner$.sketch(Partitioner.scala:309)
            at org.apache.spark.RangePartitioner.<init>(Partitioner.scala:171)
            at org.apache.spark.RangePartitioner.<init>(Partitioner.scala:151)
            at org.apache.spark.rdd.OrderedRDDFunctions$$anonfun$sortByKey$1.apply(OrderedRDDFunctions.scala:62)
            at org.apache.spark.rdd.OrderedRDDFunctions$$anonfun$sortByKey$1.apply(OrderedRDDFunctions.scala:61)
            at org.apache.spark.rdd.RDDOperationScope$.withScope(RDDOperationScope.scala:151)
            at org.apache.spark.rdd.RDDOperationScope$.withScope(RDDOperationScope.scala:112)
            at org.apache.spark.rdd.RDD.withScope(RDD.scala:385)
            at org.apache.spark.rdd.OrderedRDDFunctions.sortByKey(OrderedRDDFunctions.scala:61)
            at org.apache.spark.api.java.JavaPairRDD.sortByKey(JavaPairRDD.scala:936)
            at org.broadinstitute.hellbender.utils.spark.SparkUtils.sortUsingElementsAsKeys(SparkUtils.java:165)
            at org.broadinstitute.hellbender.utils.spark.SparkUtils.sortReadsAccordingToHeader(SparkUtils.java:143)
            at org.broadinstitute.hellbender.utils.spark.SparkUtils.querynameSortReadsIfNecessary(SparkUtils.java:306)
            at org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSpark.mark(MarkDuplicatesSpark.java:208)
            at org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSpark.mark(MarkDuplicatesSpark.java:272)
            at org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSpark.runTool(MarkDuplicatesSpark.java:354)
            at org.broadinstitute.hellbender.engine.spark.GATKSparkTool.runPipeline(GATKSparkTool.java:546)
            at org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram.doWork(SparkCommandLineProgram.java:31)
            at org.broadinstitute.hellbender.cmdline.CommandLineProgram.runTool(CommandLineProgram.java:140)
            at org.broadinstitute.hellbender.cmdline.CommandLineProgram.instanceMainPostParseArgs(CommandLineProgram.java:192)
            at org.broadinstitute.hellbender.cmdline.CommandLineProgram.instanceMain(CommandLineProgram.java:211)
            at org.broadinstitute.hellbender.Main.runCommandLineProgram(Main.java:160)
            at org.broadinstitute.hellbender.Main.mainEntry(Main.java:203)
            at org.broadinstitute.hellbender.Main.main(Main.java:289)
    25/10/15 22:41:42 INFO ShutdownHookManager: Shutdown hook called
    25/10/15 22:41:42 INFO ShutdownHookManager: Deleting directory /tmp/spark-1bc96c24-9c84-4d2e-bd64-c57bce3dfe15

    ```
</details>

## Reason 
```Unsupported class file major version 5``` means:

    - Some part of the software (in this case, Apache Spark 2.4.5, bundled inside GATK 4.3.0.0)
    was compiled with a newer Java version (Java 11),

    - But another part of your runtime environment (Spark’s internal ASM library)
    only understands older Java bytecode (Java 8).

So, Spark 2.4.5 can’t handle the .class files that come from your Java 11 environment.
It expects Java 8 (major version 52), but it’s seeing Java 11 (major version 55).

## Solution
1️⃣ Option 1 (Recommended for WSL):

Use the non-Spark version — same result, no Spark overhead.
```
gatk MarkDuplicates \
  -I SRR062634_sorted.bam \
  -O SRR062634_dedup_reads.bam \
  -M SRR062634_metrics.txt
```


It’s slower for huge datasets, but more stable and works perfectly in WSL or local machines.

Spark isn’t needed unless you’re working on a cluster or massive genomes.

2️⃣ Option 3: Use a newer GATK release with newer Spark

Newer GATK versions (4.4.x+) use Spark 3.x, which supports Java 11.

You can install that from bioconda:
```
conda install -c bioconda gatk4=4.4.0.0
```


Then re-run your Spark command — it should now work fine with Java 11.
