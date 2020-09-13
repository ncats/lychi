import os
import signal
import subprocess
import sys

def grabOutput(syscall, limit = 10000):
    pro = subprocess.Popen(syscall, stdout=subprocess.PIPE,
                           shell=True, preexec_fn=os.setsid)

    line = pro.stdout.readline()
    output = []
    while line != "" and len(output) < limit:
        output.append(line)
        if len(output)%1000 == 0:
            sys.stderr.write("lines: "+str(len(output))+"\n")
        line = pro.stdout.readline()
    try:
        os.killpg(os.getpgid(pro.pid), signal.SIGTERM)
    except:
        pass

    return output

if __name__=="__main__":

    orig = "java -cp ../lychi-comp/target/lychi-0.7.1-jar-with-dependencies.jar lychi.tautomers.SayleDelanyTautomerGenerator"
    new = '/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/bin/java "-javaagent:/Applications/IntelliJ IDEA.app/Contents/lib/idea_rt.jar=60171:/Applications/IntelliJ IDEA.app/Contents/bin" -Dfile.encoding=UTF-8 -classpath /Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/jre/lib/charsets.jar:/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/jre/lib/ext/cldrdata.jar:/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/jre/lib/ext/dnsns.jar:/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/jre/lib/ext/jaccess.jar:/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/jre/lib/ext/jfxrt.jar:/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/jre/lib/ext/localedata.jar:/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/jre/lib/ext/nashorn.jar:/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/jre/lib/ext/sunec.jar:/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/jre/lib/ext/sunjce_provider.jar:/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/jre/lib/ext/sunpkcs11.jar:/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/jre/lib/ext/zipfs.jar:/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/jre/lib/jce.jar:/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/jre/lib/jfxswt.jar:/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/jre/lib/jsse.jar:/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/jre/lib/management-agent.jar:/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/jre/lib/resources.jar:/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/jre/lib/rt.jar:/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/lib/ant-javafx.jar:/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/lib/dt.jar:/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/lib/javafx-mx.jar:/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/lib/jconsole.jar:/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/lib/packager.jar:/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/lib/sa-jdi.jar:/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/lib/tools.jar:/Users/southalln/git/lychi/binaries:/Users/southalln/git/lychi/lib/jchem3.jar:/Users/southalln/git/lychi/lib/molwitch-0.5.9.1.jar:/Users/southalln/git/lychi/lib/molwitch-jchem-20.3.0.jar:/Users/southalln/.m2/repository/gov/nih/ncats/ncats-common/0.3.1/ncats-common-0.3.1.jar lychi.tautomers.SayleDelanyTautomerGenerator'
    #new = "java -cp target/lychi-0.7.1-jar-with-dependencies.jar:lib/jchem3.jar lychi.tautomers.SayleDelanyTautomerGenerator"

    testfiles = ["ruili_all_approved_02032017_noel.smi", "9b4_Prous-pubchem_can.smi", "collisions.smi"]
    #testfiles = ['test.smiles']
    testfiles = ["ruili_all_approved_02032017_noel.smi"]

    for testfile in testfiles:
        syscall = orig + " " + testfile + " 2> /dev/null"
        output1 = grabOutput(syscall)
        print "finished run1 of "+testfile+", lines: "+str(len(output1))
        syscall = new + " " + testfile + " 2> /dev/null"
        output2 = grabOutput(syscall)
        missing = 0
        print len(output1), len(output2)
        for i in range(max(len(output1), len(output2))):
            if len(output1) > i:
                if len(output2) > i:
                    if output1[i] != output2[i]:
                        if '[S@' in output1[i] or '[S@' in output2[i]:
                            pass
                        else:
                            sline1 = output1[i].split()
                            sline2 = output2[i].split()
                            if sline1[-1] != sline2[-1]:
                                print "0", i, output1[i], output2[i]
                            else:
                                #print "1", i, output1[i], output2[i] TODO look at why the smiles output is different
                                pass
                else:
                    #print "2", i, output1[i]
                    missing = missing + 1
            else:
                #print "3", i, output2[i]
                missing = missing + 1
        print "missing", missing
