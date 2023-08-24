OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
rz(pi/4) q[1];
rz(-0.18374584309244257) q[2];
sx q[2];
rz(-1.702767829786925) q[2];
sx q[2];
rz(2.7061073589921776) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
rz(pi/2) q[4];
sx q[4];
rz(pi/2) q[4];
cx q[0],q[6];
cx q[6],q[0];
cx q[0],q[6];
rz(-pi/2) q[0];
cx q[1],q[0];
rz(-pi/2) q[0];
sx q[0];
rz(0.5047225491632941) q[0];
sx q[0];
rz(3*pi/4) q[7];
sx q[7];
rz(3*pi/4) q[7];
cx q[3],q[7];
rz(pi/4) q[7];
sx q[7];
cx q[1],q[7];
rz(-1.4426445078456234) q[1];
sx q[1];
rz(-0.4886866698351717) q[1];
rz(-pi) q[7];
sx q[7];
rz(-pi) q[7];
rz(pi/2) q[8];
cx q[5],q[8];
sx q[8];
rz(1.8636664803362741) q[8];
sx q[8];
rz(-pi) q[8];
cx q[5],q[8];
rz(1.482616996508109) q[5];
sx q[5];
rz(-1.2400174276280929) q[5];
sx q[5];
rz(-2.2886384491343073) q[5];
sx q[8];
rz(-1.8636664803362741) q[8];
sx q[8];
rz(pi/2) q[8];
cx q[6],q[8];
rz(0.846224050758152) q[8];
cx q[6],q[8];
rz(-pi/4) q[6];
rz(-pi) q[8];
sx q[8];
rz(1.4987726727285118) q[8];
sx q[8];
cx q[3],q[8];
sx q[8];
rz(1.4987726727285118) q[8];
sx q[8];
rz(-pi) q[8];
cx q[3],q[8];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
rz(3*pi/4) q[9];
sx q[9];
rz(pi/2) q[9];
cx q[4],q[9];
cx q[2],q[4];
rz(-2.5957795031712) q[4];
cx q[2],q[4];
rz(0.457168204761067) q[2];
sx q[2];
rz(-1.4358789277334028) q[2];
sx q[2];
rz(2.625442821665123) q[2];
rz(2.595779503171199) q[4];
cx q[4],q[5];
sx q[5];
rz(1.600388720916829) q[5];
sx q[5];
rz(-pi) q[5];
cx q[4],q[5];
cx q[8],q[5];
rz(0.0626643963436806) q[5];
cx q[8],q[5];
cx q[3],q[8];
rz(-pi) q[5];
sx q[5];
rz(3*pi/4) q[5];
cx q[8],q[3];
cx q[3],q[8];
rz(pi/2) q[9];
sx q[9];
rz(pi/2) q[9];
cx q[9],q[0];
sx q[0];
rz(0.5047225491632932) q[0];
sx q[0];
rz(-pi) q[0];
cx q[9],q[0];
cx q[0],q[7];
rz(2.76981809613559) q[7];
cx q[0],q[7];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4505688580845106) q[7];
rz(pi/2) q[9];
sx q[9];
rz(pi/2) q[9];
cx q[6],q[9];
rz(-pi/4) q[9];
cx q[6],q[9];
rz(0.8827324613923428) q[6];
cx q[6],q[4];
rz(-0.8827324613923437) q[4];
cx q[6],q[4];
rz(-0.6880638654025537) q[4];
sx q[4];
rz(pi/2) q[4];
cx q[4],q[1];
sx q[1];
rz(1.929877870828843) q[1];
sx q[1];
rz(-pi) q[1];
sx q[4];
rz(1.929877870828843) q[4];
sx q[4];
rz(-pi) q[4];
cx q[4],q[1];
rz(-1.0820325535127813) q[1];
sx q[1];
rz(3*pi/4) q[1];
cx q[2],q[1];
rz(-pi/4) q[1];
cx q[3],q[1];
rz(pi/4) q[1];
cx q[2],q[1];
rz(pi/4) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(3*pi/4) q[3];
sx q[3];
rz(-3*pi/4) q[3];
rz(-pi) q[4];
sx q[4];
rz(-pi) q[4];
rz(pi/4) q[6];
cx q[6],q[0];
rz(-pi/4) q[0];
cx q[6],q[0];
rz(3*pi/4) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[6];
cx q[4],q[6];
cx q[6],q[4];
rz(pi/2) q[4];
sx q[4];
rz(pi/2) q[4];
rz(pi/2) q[6];
sx q[6];
rz(pi/2) q[6];
rz(3*pi/4) q[9];
sx q[9];
rz(-pi) q[9];
cx q[9],q[5];
rz(-pi/4) q[5];
sx q[5];
cx q[5],q[7];
rz(-pi) q[5];
sx q[5];
rz(1.862208244020401) q[5];
sx q[5];
sx q[7];
rz(1.862208244020401) q[7];
sx q[7];
rz(-pi) q[7];
cx q[5],q[7];
x q[5];
rz(pi/4) q[5];
cx q[1],q[5];
rz(-pi/4) q[5];
rz(0.12022746871038636) q[7];
rz(pi/2) q[9];
sx q[9];
rz(pi/2) q[9];
cx q[0],q[9];
rz(-pi/4) q[9];
cx q[8],q[9];
rz(pi/4) q[9];
cx q[0],q[9];
rz(pi/4) q[0];
rz(-pi/4) q[9];
cx q[8],q[9];
cx q[8],q[0];
rz(-pi/4) q[0];
rz(pi/4) q[8];
cx q[8],q[0];
rz(2.854712468909513) q[0];
cx q[0],q[6];
rz(-2.8547124689095127) q[6];
cx q[0],q[6];
rz(-pi) q[0];
sx q[0];
rz(2.2656990071838656) q[0];
sx q[0];
cx q[2],q[0];
sx q[0];
rz(2.2656990071838647) q[0];
sx q[0];
rz(-pi) q[0];
cx q[2],q[0];
rz(1.0423664421441003) q[0];
sx q[0];
rz(-0.021964937464566958) q[0];
sx q[0];
rz(1.5813738674224833) q[0];
rz(pi/2) q[2];
sx q[2];
rz(3.70254445009021) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-1.8576765114751774) q[6];
sx q[6];
rz(pi/2) q[6];
cx q[8],q[5];
rz(pi/4) q[5];
cx q[1],q[5];
rz(-0.2402569818940492) q[1];
rz(pi/4) q[5];
sx q[5];
rz(pi/2) q[5];
rz(pi/4) q[8];
cx q[8],q[5];
rz(-pi/4) q[5];
cx q[8],q[5];
rz(pi/4) q[5];
cx q[5],q[1];
sx q[1];
rz(0.9869719357034512) q[1];
sx q[1];
rz(-pi) q[1];
cx q[5],q[1];
sx q[1];
rz(-0.9869719357034512) q[1];
sx q[1];
rz(-0.577761311847091) q[1];
rz(1.1947397098635086) q[8];
rz(3*pi/4) q[9];
sx q[9];
rz(pi/2) q[9];
cx q[9],q[4];
rz(-pi/4) q[4];
cx q[7],q[4];
rz(pi/4) q[4];
cx q[9],q[4];
rz(-pi/4) q[4];
cx q[7],q[4];
rz(3*pi/4) q[4];
sx q[4];
rz(-0.7112438967769421) q[4];
rz(pi/4) q[9];
cx q[7],q[9];
rz(pi/4) q[7];
rz(-pi/4) q[9];
cx q[7],q[9];
rz(1.7974487451472427) q[7];
cx q[7],q[4];
rz(0.37149643735971827) q[4];
sx q[4];
rz(-2.6624280714241655) q[4];
sx q[4];
cx q[7],q[4];
rz(1.527101694642056) q[4];
sx q[4];
rz(-2.78390337524257) q[4];
sx q[4];
rz(-2.135529099174553) q[4];
rz(0.588953893180669) q[7];
cx q[8],q[7];
rz(1.9468529437262845) q[7];
sx q[7];
rz(-0.47398457121814985) q[7];
sx q[7];
cx q[8],q[7];
rz(-2.123117044087052) q[7];
sx q[7];
rz(-1.3078797954018775) q[7];
sx q[7];
rz(1.1844745712583098) q[7];
cx q[8],q[7];
rz(3*pi/4) q[7];
sx q[7];
rz(-3*pi/4) q[7];
cx q[0],q[7];
rz(-pi/4) q[7];
rz(pi/2) q[9];
sx q[9];
rz(pi/2) q[9];
cx q[9],q[6];
rz(5.15424161688887) q[6];
cx q[9],q[6];
rz(pi/2) q[6];
sx q[6];
rz(pi/2) q[6];
cx q[3],q[6];
rz(-pi/4) q[6];
cx q[3],q[6];
cx q[2],q[3];
rz(pi/2) q[6];
cx q[6],q[5];
rz(-pi/4) q[5];
cx q[6],q[5];
rz(pi/4) q[5];
cx q[5],q[7];
rz(-pi) q[6];
sx q[6];
rz(pi/2) q[6];
cx q[4],q[6];
rz(2.47218286470388) q[6];
cx q[4],q[6];
rz(-pi) q[4];
sx q[4];
rz(-pi) q[4];
rz(-pi) q[6];
sx q[6];
rz(-pi/2) q[6];
rz(pi/4) q[7];
cx q[0],q[7];
rz(pi/4) q[7];
sx q[7];
rz(2.909385099224046) q[7];
cx q[4],q[7];
rz(-0.553190609031703) q[7];
cx q[4],q[7];
rz(2.016953725302102) q[4];
cx q[8],q[3];
rz(-pi/4) q[3];
cx q[2],q[3];
rz(pi/4) q[3];
cx q[8],q[3];
rz(-pi/4) q[3];
cx q[2],q[3];
rz(pi/4) q[3];
rz(pi/4) q[8];
cx q[2],q[8];
rz(pi/4) q[2];
rz(-pi/4) q[8];
cx q[2],q[8];
rz(-pi/2) q[2];
cx q[3],q[2];
rz(-2.099124774493133) q[2];
sx q[2];
rz(-0.5374914680525578) q[2];
sx q[2];
rz(2.775109655820465) q[2];
cx q[3],q[6];
sx q[6];
rz(0.7507350708951996) q[6];
sx q[6];
rz(-pi) q[6];
cx q[3],q[6];
rz(3.43648187643247) q[3];
sx q[3];
rz(8.90367263978313) q[3];
sx q[3];
rz(12.2526271558776) q[3];
sx q[6];
rz(-1.5361332342926488) q[6];
rz(pi/2) q[8];
sx q[8];
rz(pi/4) q[8];
rz(-pi/2) q[9];
sx q[9];
rz(-2.3814683227721294) q[9];
cx q[9],q[1];
rz(-0.6901862759156994) q[1];
sx q[1];
rz(-2.36101004398943) q[1];
sx q[1];
cx q[9],q[1];
rz(1.3523351305049989) q[1];
sx q[1];
rz(-1.3593544297422575) q[1];
sx q[1];
rz(-1.5988858292782062) q[1];
cx q[0],q[1];
rz(-pi/4) q[1];
cx q[5],q[1];
rz(pi/4) q[1];
cx q[0],q[1];
rz(pi/4) q[1];
sx q[1];
rz(pi/2) q[1];
reset q[1];
rz(1.90397825362098) q[1];
rz(pi/2) q[5];
sx q[5];
rz(pi/2) q[5];
cx q[0],q[5];
rz(-pi/4) q[5];
cx q[0],q[5];
x q[0];
cx q[0],q[5];
rz(-pi/4) q[5];
cx q[0],q[5];
cx q[0],q[6];
rz(pi/2) q[5];
sx q[5];
rz(pi/2) q[5];
cx q[4],q[5];
rz(-2.016953725302102) q[5];
cx q[4],q[5];
rz(2.8023518886995493) q[5];
cx q[5],q[4];
rz(-pi/4) q[4];
cx q[5],q[4];
rz(pi/4) q[4];
rz(pi/4) q[6];
sx q[6];
rz(pi/2) q[6];
cx q[7],q[6];
rz(pi/4) q[6];
cx q[8],q[6];
rz(-pi/4) q[6];
cx q[7],q[6];
cx q[3],q[7];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
rz(pi/4) q[6];
cx q[7],q[3];
rz(-pi/4) q[3];
cx q[8],q[6];
rz(pi/4) q[6];
sx q[6];
rz(3*pi/4) q[6];
cx q[0],q[6];
rz(-1.3387530830077825) q[0];
cx q[0],q[2];
rz(-1.8028395705820106) q[2];
sx q[2];
rz(-0.6869513817700614) q[2];
sx q[2];
cx q[0],q[2];
rz(1.5398517284786077) q[2];
sx q[2];
rz(-1.5961683227711339) q[2];
sx q[2];
rz(0.09883939712484491) q[2];
cx q[0],q[2];
rz(-pi/4) q[2];
cx q[4],q[2];
rz(pi/4) q[2];
cx q[0],q[2];
rz(3*pi/4) q[2];
sx q[2];
rz(2.054053556959518) q[4];
rz(pi/4) q[6];
sx q[6];
rz(-2.1699664750954373) q[6];
cx q[8],q[3];
rz(pi/4) q[3];
cx q[7],q[3];
rz(-pi/4) q[3];
rz(pi/4) q[7];
cx q[8],q[3];
rz(3*pi/4) q[3];
sx q[3];
rz(pi/2) q[3];
cx q[8],q[7];
rz(-pi/4) q[7];
rz(pi/4) q[8];
cx q[8],q[7];
cx q[3],q[7];
rz(-2.7083187458828704) q[3];
rz(2.45642696224234) q[7];
cx q[5],q[7];
rz(-2.45642696224234) q[7];
cx q[5],q[7];
sx q[5];
rz(2.8288941224398823) q[5];
sx q[5];
rz(-0.4420858533854628) q[7];
sx q[7];
rz(-1.4782175629780872) q[7];
sx q[7];
rz(2.4812825891826726) q[7];
rz(pi/2) q[8];
sx q[8];
rz(pi/2) q[8];
cx q[9],q[1];
rz(-1.90397825362098) q[1];
cx q[9],q[1];
cx q[6],q[1];
rz(-2.5424225052892537) q[1];
cx q[6],q[1];
rz(-2.169966475095437) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[6];
sx q[6];
rz(-0.1496276378197594) q[6];
cx q[5],q[6];
rz(0.681472444907383) q[6];
cx q[5],q[6];
rz(0.9640167043681815) q[5];
sx q[5];
rz(pi/2) q[5];
rz(pi/2) q[6];
sx q[6];
rz(pi/2) q[6];
cx q[6],q[7];
cx q[7],q[6];
cx q[6],q[7];
rz(pi/2) q[6];
cx q[8],q[1];
rz(5.9447586947432) q[1];
cx q[8],q[1];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
cx q[0],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
cx q[1],q[0];
rz(-pi/4) q[0];
cx q[3],q[0];
rz(pi/4) q[0];
cx q[1],q[0];
rz(-pi/4) q[0];
rz(pi/4) q[1];
cx q[3],q[0];
rz(3*pi/4) q[0];
sx q[0];
rz(pi/2) q[0];
cx q[3],q[1];
rz(-pi/4) q[1];
rz(pi/4) q[3];
cx q[3],q[1];
cx q[0],q[1];
sx q[1];
rz(pi/4) q[3];
cx q[3],q[0];
rz(-pi/4) q[0];
cx q[3],q[0];
rz(-2.2665017871669715) q[0];
cx q[5],q[0];
rz(-1.66048902982027) q[0];
cx q[5],q[0];
rz(pi/2) q[8];
sx q[8];
rz(-pi) q[8];
cx q[2],q[8];
cx q[8],q[2];
rz(0.9981616612409976) q[2];
sx q[2];
rz(-1.8157031069065201) q[2];
sx q[2];
rz(1.880795518825602) q[2];
cx q[3],q[8];
rz(-pi/4) q[8];
cx q[1],q[8];
rz(pi/4) q[8];
cx q[3],q[8];
rz(pi/4) q[3];
rz(-pi/4) q[8];
cx q[1],q[8];
cx q[1],q[3];
rz(pi/4) q[1];
rz(-pi/4) q[3];
cx q[1],q[3];
rz(-pi/4) q[1];
cx q[1],q[7];
reset q[3];
rz(1.3032813562413388) q[3];
rz(pi/4) q[7];
cx q[1],q[7];
rz(3*pi/4) q[1];
sx q[1];
rz(3*pi/4) q[1];
rz(pi/4) q[7];
sx q[7];
rz(pi/2) q[7];
rz(-3*pi/4) q[8];
sx q[8];
rz(pi/2) q[8];
cx q[2],q[8];
rz(0.262648434448805) q[8];
cx q[2],q[8];
x q[2];
rz(-pi/4) q[2];
rz(pi/2) q[8];
sx q[8];
cx q[7],q[8];
rz(1.47529214971582) q[8];
cx q[7],q[8];
sx q[8];
rz(-1.156895125489605) q[8];
sx q[8];
cx q[5],q[8];
sx q[8];
rz(0.4139012013052912) q[8];
sx q[8];
rz(-pi) q[8];
cx q[5],q[8];
sx q[5];
rz(0.968519349396268) q[9];
cx q[4],q[9];
rz(-0.032870309931802844) q[9];
sx q[9];
rz(-0.1073545434645844) q[9];
sx q[9];
cx q[4],q[9];
sx q[9];
rz(-0.1073545434645844) q[9];
sx q[9];
rz(2.9913417775227753) q[9];
cx q[9],q[4];
rz(-pi/4) q[4];
cx q[9],q[4];
rz(2.439777608488563) q[4];
sx q[4];
rz(pi/2) q[4];
cx q[4],q[6];
cx q[6],q[4];
cx q[0],q[4];
cx q[4],q[0];
cx q[0],q[4];
x q[0];
rz(pi/2) q[6];
sx q[6];
rz(pi/2) q[6];
cx q[3],q[6];
rz(-1.3032813562413397) q[6];
cx q[3],q[6];
cx q[3],q[7];
rz(1.3032813562413388) q[6];
cx q[6],q[1];
rz(pi/4) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/4) q[7];
reset q[9];
rz(-pi) q[9];
sx q[9];
rz(-pi) q[9];
cx q[9],q[2];
rz(pi/4) q[2];
sx q[2];
cx q[2],q[1];
rz(pi/4) q[1];
cx q[4],q[1];
rz(-pi/4) q[1];
cx q[2],q[1];
rz(pi/4) q[1];
reset q[2];
cx q[4],q[1];
rz(pi/4) q[1];
sx q[1];
rz(3*pi/4) q[1];
rz(pi/2) q[4];
sx q[4];
rz(pi/2) q[4];
cx q[6],q[1];
rz(3*pi/4) q[1];
sx q[1];
rz(pi/2) q[1];
cx q[1],q[5];
rz(5.77533654825107) q[5];
cx q[1],q[5];
x q[1];
rz(-1.5265147469991232) q[5];
sx q[5];
rz(-pi) q[5];
rz(2.447995759561678) q[6];
cx q[9],q[7];
rz(pi/4) q[7];
cx q[3],q[7];
rz(pi/4) q[3];
rz(-pi/4) q[7];
cx q[9],q[7];
rz(3*pi/4) q[7];
sx q[7];
cx q[8],q[7];
rz(pi/2) q[7];
cx q[6],q[7];
rz(pi/2) q[6];
sx q[6];
rz(pi/2) q[6];
cx q[7],q[6];
rz(-pi/4) q[6];
cx q[2],q[6];
rz(pi/4) q[6];
cx q[7],q[6];
rz(-pi/4) q[6];
cx q[2],q[6];
rz(3*pi/4) q[6];
sx q[6];
rz(pi/2) q[6];
rz(pi/4) q[7];
cx q[2],q[7];
rz(pi/4) q[2];
rz(-pi/4) q[7];
cx q[2],q[7];
cx q[6],q[7];
rz(pi/4) q[8];
cx q[9],q[3];
rz(-pi/4) q[3];
rz(pi/4) q[9];
cx q[9],q[3];
cx q[0],q[3];
rz(2.79335211007819) q[3];
cx q[0],q[3];
cx q[0],q[5];
sx q[5];
rz(1.5265147469991227) q[5];
sx q[5];
rz(-pi) q[5];
cx q[0],q[5];
rz(pi/2) q[9];
sx q[9];
rz(pi/2) q[9];
cx q[9],q[4];
rz(3.93443991868747) q[4];
cx q[9],q[4];
rz(pi/2) q[4];
sx q[4];
rz(-2.565744826837525) q[4];
rz(pi/2) q[9];
sx q[9];
rz(pi/2) q[9];
cx q[3],q[9];