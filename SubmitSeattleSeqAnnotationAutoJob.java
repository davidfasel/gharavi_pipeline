// SubmitSeattleSeqAnnotationAutoJob.java
// Peggy Robertson, Nickerson Lab, Dept. of Genome Sciences, University of Washington
// 1/2014
// example client program for querying SeattleSeq Annotation from a local computer; runs with Java 1.5 or higher

/*

steps for using:

1. A line like this must be in the file submitted; the filename is only used by the server, and does not need to be unique (a time-stamp is added by the server):
# autoFile testAuto.txt

2. If you want to receive gzipped (compressed) files from the server, add this line to the input file (recommended for large jobs):
# compressAuto true

3. Modify this code to put your own inputFile name, outputFile name, eMail, compression choice, and timeout minutes in.    In the submitTheInputFile function, choose annotation parameters.
     The code will not compile until an email address is added.

4. Execute "java -version" to make sure java is installed and has a version of at least 1.5 (plus any subversion).

5. Acquire some jar files if not present (newer versions may be ok):
    httpunit.jar
    nekohtml-0.9.5 or jtidy-4aug2000r7-dev.jar
    js-1.6R5.jar
    xercesImpl-2.6.1.jar

6. Set your classpath variable with your own path substituted for "path" below (this is for bash, and can be put in your .bashrc file so you only have to do it once):
export CLASSPATH=./:/path/httpunit.jar:/path/js-1.6R5.jar:/path/nekohtml-0.9.5.jar:/path/xercesImpl.jar

7. Compile:
javac SubmitSeattleSeqAnnotationAutoJob.java

8. Execute:
java SubmitSeattleSeqAnnotationAutoJob

If your job fails and you see a line like "downloadURL = <html>", you've probably forgotten to put the autoFile line in your submitted file.

*/

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.BufferedWriter;
import java.io.PrintWriter;
import java.io.FileWriter;
import java.io.FileOutputStream;

import java.util.zip.GZIPInputStream;

//import junit.framework.*;
import com.meterware.httpunit.*;
import com.meterware.httpunit.protocol.*;

//public class SubmitSeattleSeqAnnotationAutoJob extends TestCase {
public class SubmitSeattleSeqAnnotationAutoJob {
        // ----------------------------- the parameters in this section must be customized
    int timeoutMinutes = 1000;    // increment this for long jobs
    private String eMail = "daf2139@columbia.edu";    // must be filled in; an email will not normally be sent, but error messages will be
    private boolean useCompressionForOutput = true;    // add a compressAuto line in the input file if true

    private String inputFilename, outputFilename;
    final int BUFFER_SIZE = 65536;
    private String URL = "http://snp.gs.washington.edu/SeattleSeqAnnotation138/";    // SeattleSeq Annotation URL
    private String downloadURL = "";    // URL to get result file
    private String progressURL = "";    // URL to monitor when job is done
    WebConversation wc = null;

    public SubmitSeattleSeqAnnotationAutoJob(String argv[]) {
            if (argv.length < 2) {
                System.err.println("Please specify an input and output file path.    Usage:");
                System.err.println("java SubmitSeattleSeqAnnotationAutoJob inputFilePath outputFilePath");
                System.exit(0);
        }
        inputFilename = argv[0];
        outputFilename = argv[1];
    }

    public void runJob() {
        try {
            // create the conversation object
            wc = new WebConversation();

             // obtain the main page on the SeattleSeq Annotation site
            System.out.println("     Asking for the main SeaSeq page...");
            WebResponse response = getMainPageResponse(wc);

            // submit the request
            System.out.print("     Submitting file " + inputFilename + "...");
            // fill out the form on the main page with the file, the email address, and the options
            response = submitTheInputFile(response, inputFilename, wc);

            // set the URLs for download and progress
            boolean validURLs = setURLsFromSource(response);
            if (!validURLs) {
                System.out.println("The attempt to get URLs for monitoring progress failed.    " +
                    "Most likely your submitted file does not have an autoFile line in it, or is not of the expected format.");
            }
            else if (downloadURL.contains("ABORT") || progressURL.contains("ABORT")) {
                System.out.println("ERROR: The SeattleSeqAnnotation file was aborted.    \n" +
                      "The job may be too large or too many simultaneous jobs may have been submitted." +
                      "Download URL:\n" + downloadURL + "\n" +
                      "Progress URL:\n " + progressURL + "\n" +
                      "To cancel job use this link:\n" + 
                      "http://snp.gs.washington.edu/SeattleSeqAnnotation137/BatchCancelServlet?cancelFile=AnnotationCancel. + jobid");
            }
            else {    // the normal case
                // e.g.
                // http://gvsbatch.gs.washington.edu/SeattleSeqAnnotation138/BatchFileDownloadServlet?file=testAuto.123456789.txt&download=plainText
                // http://gvsbatch.gs.washington.edu/SeattleSeqAnnotation138/BatchProgressServlet?progressFile=progress.testAuto.123456789.txt&auto=true
                System.out.println(" Success! \n\n" +
                        "Download URL:\n " + downloadURL + "\n" +
                        "Progress URL:\n " + progressURL + "\n");

                // wait for the job to finish, monitoring with the progress file
                boolean isSuccessfulFinish = waitForJobToFinish();    // use progressURL to detect job completion

                // download the result file
                if (isSuccessfulFinish) {
                    WebRequest requestDownload = new GetMethodWebRequest(downloadURL);    // SeattleSeq Annotation uses http GET here
                    WebResponse responseDownload = wc.getResponse(requestDownload);
                    writeSourceToFile(responseDownload);    // write the downloaded results to a local file
                    System.out.println("download complete");
                }
                else {
                    System.out.println("PROBLEM: could not detect successful finish of job");
                }
            }
        } catch (Exception ex) {
            System.err.println("ERROR runJob - " + ex);
            ex.printStackTrace();
        }
    }

    private WebResponse getMainPageResponse(WebConversation wc) {
        // obtain the main page on the SeattleSeq Annotation site
        WebResponse responseResult = null;
        try {
            WebRequest request = new GetMethodWebRequest(URL);
            responseResult = wc.getResponse(request);
        } catch (Exception ex) {
            System.err.println("ERROR getMainPageResponse - " + ex);
        }
        return responseResult;
    }

    private WebResponse submitTheInputFile(WebResponse response, String inputFilename, WebConversation wc) {
        WebResponse responseInput = null;
        try {
            // fill out the form for submitting a file
            WebForm form = response.getFormWithName("GenotypeSummary");                // GenotypeSummary is the name of the form on SeattleSeq Annotation

            String contentType = "";
            int i = inputFilename.lastIndexOf('.');
            if (i > 0) {
                    String extension = inputFilename.substring(i+1);
                    if (extension.equals("gz")) {contentType = "application/x-gzip";}
            }

            UploadFileSpec fileSpec = new UploadFileSpec(new File(inputFilename), contentType);
            form.setParameter("GenotypeFile", new UploadFileSpec[] { fileSpec });    // GenotypeFile is the name of the file field
            form.setParameter("EMail", eMail);                                        // EMail is the name of the email field
            //form.getScriptableObject().setParameterValue( "autoFile", "UW" );

            // set various options; more parameter names may be found by viewing the html source for the home page

            /* Input File Format */
            //form.setParameter("fileFormat", "Maq");    // use this to choose an input format that is not VCF SNPs
            form.setParameter("outputFileFormatBoth", "VCFOutBoth");

            /* Additional Annotations */

            // un-comment the next two lines if your individual has a dbSNP ID, and change 5145 to the appropriate value
            // this must be set before "individual" is set (if used)
            //form.setCheckbox("columns", "dbSNPGenotype", true);
            //form.setParameter("individual", "5145");

            // CADD score: have to specifically turn this on, it's off by default
            form.setCheckbox("columns", "scoreCADD", true);

            /* DAF: turning a few of these off that aren't especially helpful in
                 hopes of speeding up the annotation a bit. */

            // Alleles in dbSNP (just gives a true/false value)
            form.setCheckbox("columns", "allelesDBSNP", false);

            // HapMap Frequencies (very few were annotated and we have frequencies from Annovar, etc)
            form.setCheckbox("columns", "HapMapFreq", false);

            // Has Genotypes
            form.setCheckbox("columns", "hasGenotypes", false);

            //dbSNP Validation
            form.setCheckbox("columns", "dbSNPValidation", false);

            // turn off "NHLBI ESP Allele Counts"
            form.setCheckbox("columns", "genomesESP", false);



            // use this to choose HapMap frequencies by reference-allele
            // form.setParameter("HapMapFreqType", "HapMapFreqRef");

            form.setCheckbox("ESPType", "ESPByPop", true);

            // click the submit button
            SubmitButton button = form.getSubmitButton("gFetch");                    // gFetch is the name of the submit button
            WebRequest requestForm = form.getRequest(button);
            responseInput = wc.getResponse(requestForm);
        } catch (Exception ex) {
            System.err.println("ERROR submitTheInputFile - " + ex);
        }

        return responseInput;
    }

    private void writeSourceToFile(WebResponse response) {
        if (useCompressionForOutput) {
            System.out.println("The submitted file is expected to have a compressAuto-true line, and the returned file is expected to be gzipped.");
            writeCompressedSourceToFile(response);
        }
        else {
            System.out.println("The submitted file is expected to have no compressAuto-true line, and the returned file is expected to be plain text.");
            writeUncompressedSourceToFile(response);
        }
    }

    private void writeUncompressedSourceToFile(WebResponse response) {
        // write the results to a local file
        try {
            PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(outputFilename), BUFFER_SIZE));    // writes to local file
            InputStream inputStream = response.getInputStream();
            InputStreamReader reader = new InputStreamReader(inputStream);
            BufferedReader bufferedReader = new BufferedReader(reader);    // reads the plain-text source from the server
            String line1 = null;
            while((line1 = bufferedReader.readLine()) != null) {    // read a line, then write a line
                writer.println(line1);
            }
            inputStream.close();
            bufferedReader.close();
            writer.flush();
            writer.close();
        } catch (Exception ex) {
            System.err.println("ERROR writeUncompressedSourceToFile - " + ex);
            ex.printStackTrace();
        }
    }

    private void writeCompressedSourceToFile(WebResponse response) {
        // write the results to a local file
        try {
            String contentType = response.getContentType();
            InputStream inputStream = response.getInputStream();
            if (contentType.equals("application/x-gzip") || contentType.equals("application/x-download")) {
                String compressedOutputFilename = outputFilename;
                if (!compressedOutputFilename.endsWith(".gz")) {
                    compressedOutputFilename += ".gz";
                }
                File outputFile = new File(compressedOutputFilename);
                FileOutputStream outputStream = new FileOutputStream(outputFile);
                byte[] buf = new byte[4 * 1024];    // 4K char buffer
                int bytesRead;
                while ((bytesRead = inputStream.read(buf)) != -1) {
                    outputStream.write(buf, 0, bytesRead);
                }
                outputStream.flush();
                outputStream.close();
            }
            else {
                System.err.println("ERROR writeCompressedSourceToFile: a compressed content type for the returned file was not detected; check the # compress line of the input file");
            }
        } catch (Exception ex) {
            System.err.println("ERROR writeCompressedSourceToFile - " + ex);
            ex.printStackTrace();
        }
    }

    private boolean setURLsFromSource(WebResponse response) {
        // process the text from the submission-acknowledge http message, and extract the links for monitoring progress and downloading the result file
        boolean validResponse = true;
        try {
            InputStream inputStream = response.getInputStream();
            InputStreamReader reader = new InputStreamReader(inputStream);
            BufferedReader bufferedReader = new BufferedReader(reader);
            String line1 = bufferedReader.readLine();

            /* DAF debug: print entire response */
//             String sCurrentLine = line1;
//             System.err.println(sCurrentLine);
//             while ((sCurrentLine = bufferedReader.readLine()) != null) {
//                    System.err.println(sCurrentLine);
//             }

            // the first line in the source contains the download and progress URLs, comma-separated
            String[] parts = line1.split(",");
            // if we got HTML back then something went wrong
            if (line1.toLowerCase().contains("html")) {
                validResponse = false;
            }
            else {    // the normal case
                downloadURL = parts[0];
                progressURL = parts[1];
                if (downloadURL.contains("localhost")) {
                    downloadURL = downloadURL.replace("localhost", "doshi.gs.washington.edu");
                }
                if (progressURL.contains("localhost")) {
                    progressURL = progressURL.replace("localhost", "doshi.gs.washington.edu");
                }
            }
            inputStream.close();
            bufferedReader.close();
        } catch (Exception ex) {
            System.err.println("ERROR setURLsFromSource - " + ex);
            ex.printStackTrace();
        }
        return validResponse;
    }

    private boolean waitForJobToFinish() {
        // keep reading the progress file until the job is done, or it times out
        // the web response is two comma-separated numbers; the first is the number of lines completed, the second is the total number of lines
        // at the beginning, before any lines are processed, the response is 0,0

        boolean isSuccessfulFinish = false;
        int sleepMilliseconds = 10000;    // 10 seconds; this could be increased for long jobs
        int maxCycles = timeoutMinutes * 60000 / sleepMilliseconds;

        try {
            System.out.println("Checking progress every 10 seconds:");
            int numberTries = 0;

            while (true) {
                numberTries++;
                if (numberTries > maxCycles) {    // failed to finish in alloted time
                    System.out.println("ERROR: the annotation timed out after " + numberTries + " cycles, time-out minutes is " + timeoutMinutes);
                    break;
                }

                // read information from the progress file on the server
                // the progress URL simply returns two numbers "35700,597233" (lines finished / total lines)
                WebRequest requestProgress = new GetMethodWebRequest(progressURL);    //BatchProgressServlet uses GET
                WebResponse responseProgress = wc.getResponse(requestProgress);
                InputStreamReader reader = new InputStreamReader(responseProgress.getInputStream());
                BufferedReader bufferedReader = new BufferedReader(reader);
                String line1 = bufferedReader.readLine();
                String[] parts = line1.split(",");
                int numberLinesFinished = Integer.valueOf(parts[0]).intValue();
                int numberLinesTotal = Integer.valueOf(parts[1]).intValue();
                bufferedReader.close();
                System.out.print("Cycle:" + numberTries + ", Lines: " + numberLinesFinished + " out of total " + numberLinesTotal + "\r");

                // see if finished
                if (numberLinesFinished == numberLinesTotal && numberLinesFinished != 0) {
                    isSuccessfulFinish = true;
                    break;
                }
                Thread.sleep(sleepMilliseconds);
            }
            Thread.sleep(5000);    // give it 5 more seconds, to be sure
            System.out.println("end wait");
        } catch (Exception ex) {
            System.err.println("ERROR waitForJobToFinish - " + ex);
        }
        return isSuccessfulFinish;
    }

    public static void main(String argv[]) throws Exception {
        SubmitSeattleSeqAnnotationAutoJob myself = new SubmitSeattleSeqAnnotationAutoJob(argv);
        myself.runJob();
    }

}
