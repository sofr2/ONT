package uk.ac.ebi.ena.ont;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileSpan;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import com.beust.jcommander.converters.FileConverter;

public class Identity {
	private byte[] ref;
	private int[] idents;
	private int[] nonIdents;
	private int[] inserts;
	private int[] dels;
	private int[] depth;
	private long unmappedBaseCount;
	private long totalBaseCount;

	private long unalignedBaseCount;
	private long alignedBaseCount;

	private long identicalBaseCount;
	private long nonIdenticalBaseCount;

	public Identity(byte[] ref) {
		this.ref = ref;
		idents = new int[ref.length];
		nonIdents = new int[ref.length];
		inserts = new int[ref.length];
		dels = new int[ref.length];
		depth = new int[ref.length];
	}

	private static void printUsage(JCommander jc) {
		StringBuilder sb = new StringBuilder();
		sb.append("\n");
		jc.usage(sb);

		System.out.println("Version " + Identity.class.getPackage().getImplementationVersion());
		System.out.println(sb.toString());
	}

	public static void main(String[] args) throws IOException {

		Params params = new Params();
		JCommander jc = new JCommander(params);
		try {
			jc.parse(args);
		} catch (Exception e) {
			System.err.println("Failed to parse parameteres, detailed message below: ");
			System.err.println(e.getMessage());
			System.err.println();
			System.err.println("See usage: -h");
			System.exit(1);
		}

		if (args.length == 0 || params.help) {
			printUsage(jc);
			System.exit(1);
		}

		if (params.referenceFasta == null) {
			System.err
					.println("No reference file specified, remote access over internet may be used to download public sequences. ");
			System.exit(1);
		}

		if (params.fileList == null || params.fileList.isEmpty()) {
			System.err.println("Expecting one or more input files.");
			System.exit(1);
		}

		if (params.sequenceName == null) {
			SAMFileReader reader = new SAMFileReader(new File(params.fileList.get(0)));
			SAMRecordIterator iterator = reader.iterator();
			SAMRecord record = iterator.next();
			params.sequenceName = record.getReferenceName();
			reader.close();
		}
		ReferenceSequenceFile rsf = new IndexedFastaSequenceFile(params.referenceFasta);
		ReferenceSequence sequence;
		PrintStream ps = null;
		if (params.identityFileName != null)
			ps = "-".equals(params.identityFileName) ? ps = System.out : new PrintStream(new File(
					params.identityFileName));

		while ((sequence = rsf.nextSequence()) != null) {
			long nofRecords = 0;
			Identity identity = null;

			for (String filePath : params.fileList) {
				File bamFile = new File(filePath);
				if (!bamFile.canRead()) {
					System.err.println("Can't read file: " + bamFile.getAbsolutePath());
					System.exit(1);
				}
				SAMFileReader reader = new SAMFileReader(bamFile);
				try {

					SAMFileSpan span = reader.getIndex().getSpanOverlapping(
							reader.getFileHeader().getSequenceIndex(sequence.getName()), -1, -1);
					SAMRecordIterator iterator = reader.iterator(span);
					while (iterator.hasNext()) {
						SAMRecord record = iterator.next();

						if (params.prohibitFlags != 0 && (params.prohibitFlags & record.getFlags()) != 0)
							continue;
						if (params.requireFlags != 0
								&& (params.requireFlags & record.getFlags()) != params.requireFlags)
							continue;
						if (identity == null)
							identity = new Identity(rsf.getSequence(sequence.getName()).getBases());
						identity.addRead(record);
						nofRecords++;
					}
				} finally {
					try {
						reader.close();
					} catch (Exception e) {
					}
				}
			}

			if (identity == null)
				continue;

			float coverage;
			if (identity.totalBaseCount == 0 || identity.totalBaseCount == identity.unmappedBaseCount)
				coverage = 0;
			else
				coverage = 100f * identity.identicalBaseCount / (identity.totalBaseCount - identity.unmappedBaseCount);

			float[] series = new float[identity.inserts.length];
			for (int i = 0; i < identity.idents.length; i++) {
				series[i] = identity.idents[i] + identity.nonIdents[i] + identity.inserts[i] + identity.dels[i];
			}

			float mean = 0;
			int n = 0;
			float M2 = 0;
			float delta;
			for (int i = 0; i < series.length; i++) {
				n = n + 1;
				delta = series[i] - mean;
				mean = mean + delta / n;
				M2 = M2 + delta * (series[i] - mean);
			}

			float variance = M2 / (n - 1);
			float std = (float) Math.sqrt(variance);

			System.out.printf("@%s\t%d\t%s\t%.2f\t%.2f\t%.2f\n", sequence.getName(), nofRecords, identity.toString(),
					coverage, mean, std);

			int bins = 100;
			int[] maxBins = new int[bins];
			float[] covBins = new float[bins];
			int[] avg = new int[bins];
			for (int i = 0; i < identity.idents.length; i++) {
				int bin = (int) ((float) bins / identity.idents.length * i);
				if (maxBins[bin] < series[i])
					maxBins[bin] = (int) series[i];
			}

			for (int b = 0; b < bins; b++) {
				int binWidth = identity.idents.length / bins;
				int max = 0;
				int sum = 0;
				float noCov = 0f;
				int count = 0;

				for (int i = 0; i < binWidth; i++) {
					int p = binWidth * b + i;
					if (p >= series.length)
						break;
					count++;
					if (max < series[p])
						max = (int) series[p];
					sum += (int) series[p];
					if (series[p] < 1f)
						noCov++;
				}
				maxBins[b] = max;
				avg[b] = sum / count;
				covBins[b] = 1f - noCov / count;
			}

			for (int m : maxBins) {
				System.out.print(toChar(m, mean, std));
			}
			System.out.println();

			for (float m : covBins) {
				if (m == 0f)
					System.out.print('_');
				else {
					System.out.print((char) ('A' + m * 25f));
				}
			}
			System.out.println();

			if (ps != null) {
				for (int i = 0; i < identity.idents.length; i++) {
					int depth = identity.idents[i] + identity.nonIdents[i] + identity.inserts[i] + identity.dels[i];
					char d = (char) ('A' + depth / 1000);
					d = toChar(depth, mean, std);
					ps.println(String.format("%s\t%d\t%d\t%d\t%d\t%d\t%c", sequence.getName(), i, identity.idents[i],
							identity.nonIdents[i], identity.inserts[i], identity.dels[i], d));
				}

			}
		}
		if (ps != null)
			ps.close();
	}

	private static char toChar(int m, float mean, float std) {
		if (m >= mean) {
			int v = (int) ((m - mean) / std);
			if (v < 26)
				return (char) ('A' + (m - mean) / std);
			else
				return '^';

		} else
			return '_';
	}

	public Identity() {
		super();
	}

	@Override
	public String toString() {
		return String.format("%d\t%d\t%d\t%d\t%d\t%d",
				// @formatter:off
				unmappedBaseCount,
				totalBaseCount,
				unalignedBaseCount,
				alignedBaseCount,
				identicalBaseCount,
				nonIdenticalBaseCount
				// @formatter:on
				);
	}

	public int addRead(SAMRecord record) {
		totalBaseCount += record.getReadLength();
		if (record.getReadUnmappedFlag()) {
			unmappedBaseCount += record.getReadLength();
			return 0;
		}

		int identicalBases = 0;
		unalignedBaseCount += record.getReadLength();
		for (AlignmentBlock block : record.getAlignmentBlocks()) {
			alignedBaseCount += block.getLength();
			unalignedBaseCount -= block.getLength();
			for (int i = 0; i < block.getLength(); i++) {
				// premature cycle break ignores bases mapped beyond the
				// reference:
				if (i + block.getReferenceStart() - 1 > ref.length)
					break;

				byte base = record.getReadBases()[i + block.getReadStart() - 1];
				byte refBase = ref[i + block.getReferenceStart() - 1];
				if (base == refBase) {
					identicalBases++;
					idents[i + block.getReferenceStart() - 1]++;
				} else {
					nonIdenticalBaseCount++;
					nonIdents[i + block.getReferenceStart() - 1]++;
				}
			}
		}
		identicalBaseCount += identicalBases;

		int pos = record.getAlignmentStart();
		for (CigarElement e : record.getCigar().getCigarElements()) {
			// dumb guard against reads mapped beyond ref end:
			if (pos > inserts.length - 1)
				break;
			switch (e.getOperator()) {
			case I:
				inserts[pos]++;
				break;
			case D:
				dels[pos]++;
				break;

			default:
				break;
			}
			if (e.getOperator().consumesReferenceBases())
				pos += e.getLength();
		}

		return identicalBases;
	}

	public byte[] getRef() {
		return ref;
	}

	public long getUnmappedBaseCount() {
		return unmappedBaseCount;
	}

	public long getTotalBaseCount() {
		return totalBaseCount;
	}

	public long getUnalignedBaseCount() {
		return unalignedBaseCount;
	}

	public long getAlignedBaseCount() {
		return alignedBaseCount;
	}

	public long getIdenticalBaseCount() {
		return identicalBaseCount;
	}

	public long getNonIdenticalBaseCount() {
		return nonIdenticalBaseCount;
	}

	@Parameters(commandDescription = "Identity calulator. ")
	static class Params {
		@Parameter(names = { "--reference-fasta-file", "-R" }, converter = FileConverter.class, description = "The reference fasta file, uncompressed and indexed (.fai file, use 'samtools faidx'). ")
		File referenceFasta;

		@Parameter(names = { "--query", "-Q" })
		String sequenceName;

		@Parameter
		List<String> fileList;

		@Parameter(names = { "--help", "-h" }, arity = 0)
		boolean help = false;

		@Parameter(names = { "--identity-output", "-O" })
		String identityFileName;

		@Parameter(names = { "--prohibit-flags", "-F" })
		int prohibitFlags = 0;

		@Parameter(names = { "--require-flags", "-f" })
		int requireFlags = 0;
	}
}
