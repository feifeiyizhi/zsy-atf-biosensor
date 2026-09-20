import urllib.request, urllib.parse, ssl, time, re, json, sys
ctx=ssl.create_default_context(); ctx.check_hostname=False; ctx.verify_mode=ssl.CERT_NONE
WT="MQKKLTRSQQKHLDIINAAKEEFIEFGFLAANMDRITSSAEVSKRTLYRHFESKEVLFESVLTIINDSVNESISYHFDPNKSTEEQLTEIAYKEIDVLYKTYGIALARTIVMEFLRQPEMAKTLIQNIYSIRAITQWFRSAIEAKRLKDADPKLMTDVYVSLFQGLFFWPQVMHLDLEPHGEELSQKIETLTTIFLQSYGVAE"
params={"CMD":"Put","PROGRAM":"blastp","DATABASE":"pdb","QUERY":WT,"EXPECT":10,"HITLIST_SIZE":20,"FORMAT_TYPE":"JSON2","ALIGNMENTS":20}
data=urllib.parse.urlencode(params).encode()
for _ in range(5):
    try:
        r=urllib.request.urlopen(urllib.request.Request("https://blast.ncbi.nlm.nih.gov/Blast.cgi",data=data),timeout=60,context=ctx)
        t=r.read().decode(); m=re.search(r'RID = (\S+)',t)
        if m: rid=m.group(1); break
    except Exception as e: time.sleep(8)
if not rid: sys.exit("no RID")
open("structure/blast_rid.txt","w").write(rid)
for i in range(90):  # 最多 ~12分钟
    time.sleep(8)
    try:
        r=urllib.request.urlopen("https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_TYPE=JSON2&RID=%s"%rid,timeout=60,context=ctx)
        t=r.read().decode()
    except Exception as e: continue
    if '"status":"WAITING"' in t: continue
    if '"status":"UNKNOWN"' in t: open("structure/blast_out.json","w").write(t); sys.exit("UNKNOWN")
    if '"status":"READY"' in t or '"BlastOutput_iterations"' in t or '"hits"' in t:
        open("structure/blast_out.json","w").write(t); print("READY"); sys.exit(0)
sys.exit("timeout")
