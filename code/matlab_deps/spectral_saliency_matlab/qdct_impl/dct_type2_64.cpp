/**
 * Copyright 2011 B. Schauerte. All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions are 
 * met:
 * 
 *    1. Redistributions of source code must retain the above copyright 
 *       notice, this list of conditions and the following disclaimer.
 * 
 *    2. Redistributions in binary form must reproduce the above copyright 
 *       notice, this list of conditions and the following disclaimer in 
 *       the documentation and/or other materials provided with the 
 *       distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY B. SCHAUERTE ''AS IS'' AND ANY EXPRESS OR 
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
 * DISCLAIMED. IN NO EVENT SHALL B. SCHAUERTE OR CONTRIBUTORS BE LIABLE 
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR 
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *  
 * The views and conclusions contained in the software and documentation
 * are those of the authors and should not be interpreted as representing 
 * official policies, either expressed or implied, of B. Schauerte.
 */

/**
 * If you use any of this work in scientific research or as part of a larger
 * software system, you are kindly requested to cite the use in any related 
 * publications or technical documentation. The work is based upon:
 *
 * [1] B. Schauerte, and R. Stiefelhagen, "Predicting Human Gaze using 
 *     Quaternion DCT Image Signature Saliency and Face Detection," in IEEE 
 *     Workshop on the Applications of Computer Vision (WACV), 2012.
 * [2] B. Schauerte, and R. Stiefelhagen, "Quaternion-based Spectral 
 *     Saliency Detection for Eye Fixation Prediction," in European 
 *     Conference on Computer Vision (ECCV), 2012
 */
#include "dct_type2.hpp"

#include <cstddef>				 // for size_t

#ifndef DK
#define DK(name, value) const E name = K(value)
#endif

#ifndef WS
#define WS(stride, i)  (stride * i)
#endif

#ifndef MAKE_VOLATILE_STRIDE
#define MAKE_VOLATILE_STRIDE(x) x
#endif

template <typename R, typename stride, typename INT>
void
dct_type2_64(const R * I, R * O, stride is, stride os, INT v, INT ivs, INT ovs)
{
	typedef R E;
	typedef R K;
	DK(KP341923777, +0.341923777520602452727284714416527063932658118);
	DK(KP1_970555284, +1.970555284777882489548036866357095574320258312);
	DK(KP438202480, +0.438202480313739594455475094994715597696721593);
	DK(KP1_951404260, +1.951404260077057088920791532839055943288024532);
	DK(KP1_151616382, +1.151616382835690601491944907631461683552016911);
	DK(KP1_635169626, +1.635169626303167393009841768261267618942085035);
	DK(KP1_069995239, +1.069995239774194421326153809274035831120531384);
	DK(KP1_689707130, +1.689707130499414146519142410209914195439571963);
	DK(KP773010453, +0.773010453362736960810906609758469800971041293);
	DK(KP634393284, +0.634393284163645498215171613225493370675687095);
	DK(KP1_306345685, +1.306345685907553528168406027312610830153720047);
	DK(KP1_514417693, +1.514417693012969095150928107211568946080867431);
	DK(KP899222659, +0.899222659309213200092589158848454151766374097);
	DK(KP1_786448602, +1.786448602391030640684832894986795956001251178);
	DK(KP147129127, +0.147129127199334847058931243150468643626598531);
	DK(KP1_994580913, +1.994580913357380432271194280365135642343377358);
	DK(KP627363480, +0.627363480797782953312957691988200619986755019);
	DK(KP1_899056361, +1.899056361186073334391872148378690056504448308);
	DK(KP290284677, +0.290284677254462367636192375817395274691476278);
	DK(KP956940335, +0.956940335732208864935797886980269969482849206);
	DK(KP1_191398608, +1.191398608984866686934073057659939779023852677);
	DK(KP1_606415062, +1.606415062961289819613353025926283847759138854);
	DK(KP1_028205488, +1.028205488386443453187387677937631545216098241);
	DK(KP1_715457220, +1.715457220000544139804539968569540274084981599);
	DK(KP293460948, +0.293460948910723503317700259293435639412430633);
	DK(KP1_978353019, +1.978353019929561946903347476032486127967379067);
	DK(KP485960359, +0.485960359806527779896548324154942236641981567);
	DK(KP1_940062506, +1.940062506389087985207968414572200502913731924);
	DK(KP049082457, +0.049082457045824576063469058918565850130932238);
	DK(KP1_999397637, +1.999397637392408440231531299332344393700122163);
	DK(KP719790073, +0.719790073069976297550209144653512840404634842);
	DK(KP1_865985597, +1.865985597669477775423320511086604996590031041);
	DK(KP1_379081089, +1.379081089474133849233461259914969405691073689);
	DK(KP1_448494165, +1.448494165902933841882138486581106334966186010);
	DK(KP810482628, +0.810482628009979741816962611010104933023895508);
	DK(KP1_828419511, +1.828419511407061309270029658787154802089382231);
	DK(KP995184726, +0.995184726672196886244836953109479921575474869);
	DK(KP098017140, +0.098017140329560601994195563888641845861136673);
	DK(KP1_230463181, +1.230463181161253690969827126827968555318860016);
	DK(KP1_576692855, +1.576692855253212524018329410719378565312986274);
	DK(KP985796384, +0.985796384459568073746053377517618536479374613);
	DK(KP1_740173982, +1.740173982217422837304584808967697687821655579);
	DK(KP244821350, +0.244821350398432396997408948301891575150447218);
	DK(KP1_984959069, +1.984959069197419996313534503322235640021641309);
	DK(KP533425514, +0.533425514949796772650573030232872788084233977);
	DK(KP1_927552131, +1.927552131590879733372928711015670307326167698);
	DK(KP471396736, +0.471396736825997648556387625905254377657460319);
	DK(KP881921264, +0.881921264348355029712756863660388349508442621);
	DK(KP1_343117909, +1.343117909694036801250753700854843606457501264);
	DK(KP1_481902250, +1.481902250709918182351233794990325459457910619);
	DK(KP855110186, +0.855110186860564188641933713777597068609157259);
	DK(KP1_807978586, +1.807978586246886663172400594461074097420264050);
	DK(KP098135348, +0.098135348654836028509909953885365316629490726);
	DK(KP1_997590912, +1.997590912410344785429543209518201388886407229);
	DK(KP673779706, +0.673779706784440101378506425238295140955533559);
	DK(KP1_883088130, +1.883088130366041556825018805199004714371179592);
	DK(KP196034280, +0.196034280659121203988391127777283691722273346);
	DK(KP1_990369453, +1.990369453344393772489673906218959843150949737);
	DK(KP580569354, +0.580569354508924735272384751634790549382952557);
	DK(KP1_913880671, +1.913880671464417729871595773960539938965698411);
	DK(KP1_268786568, +1.268786568327290996430343226450986741351374190);
	DK(KP1_546020906, +1.546020906725473921621813219516939601942082586);
	DK(KP942793473, +0.942793473651995297112775251810508755314920638);
	DK(KP1_763842528, +1.763842528696710059425513727320776699016885241);
	DK(KP390180644, +0.390180644032256535696569736954044481855383236);
	DK(KP1_961570560, +1.961570560806460898252364472268478073947867462);
	DK(KP1_111140466, +1.111140466039204449485661627897065748749874382);
	DK(KP1_662939224, +1.662939224605090474157576755235811513477121624);
	DK(KP1_414213562, +1.414213562373095048801688724209698078569671875);
	DK(KP2_000000000, +2.000000000000000000000000000000000000000000000);
	DK(KP765366864, +0.765366864730179543456919968060797733522689125);
	DK(KP1_847759065, +1.847759065022573512256366378793576573644833252);
	DK(KP195090322, +0.195090322016128267848284868477022240927691618);
	DK(KP980785280, +0.980785280403230449126182236134239036973933731);
	DK(KP555570233, +0.555570233019602224742830813948532874374937191);
	DK(KP831469612, +0.831469612302545237078788377617905756738560812);
	DK(KP923879532, +0.923879532511286756128183189396788286822416626);
	DK(KP382683432, +0.382683432365089771728459984030398866761344562);
	DK(KP707106781, +0.707106781186547524400844362104849039284835938);
	INT i;
	for (i = v; i > 0; i = i - 1, I = I + ivs, O = O + ovs, MAKE_VOLATILE_STRIDE(is), MAKE_VOLATILE_STRIDE(os))
	{
		E T273;
		E T340;
		E T11;
		E T195;
		E T150;
		E T228;
		E T387;
		E T431;
		E T26;
		E T229;
		E T418;
		E T432;
		E T145;
		E T196;
		E T280;
		E T337;
		E T44;
		E T140;
		E T390;
		E T443;
		E T200;
		E T225;
		E T288;
		E T334;
		E T61;
		E T141;
		E T393;
		E T444;
		E T203;
		E T226;
		E T295;
		E T335;
		E T127;
		E T173;
		E T405;
		E T438;
		E T408;
		E T437;
		E T136;
		E T174;
		E T215;
		E T249;
		E T325;
		E T363;
		E T330;
		E T364;
		E T218;
		E T250;
		E T90;
		E T170;
		E T398;
		E T434;
		E T401;
		E T435;
		E T99;
		E T171;
		E T208;
		E T246;
		E T308;
		E T360;
		E T313;
		E T361;
		E T211;
		E T247;
		{
			E T3;
			E T271;
			E T149;
			E T272;
			E T6;
			E T338;
			E T9;
			E T339;
			{
				E T1;
				E T2;
				E T147;
				E T148;
				T1 = I[0];
				T2 = I[WS(is, 63)];
				T3 = T1 - T2;
				T271 = T1 + T2;
				T147 = I[WS(is, 32)];
				T148 = I[WS(is, 31)];
				T149 = T147 - T148;
				T272 = T147 + T148;
			}
			{
				E T4;
				E T5;
				E T7;
				E T8;
				T4 = I[WS(is, 16)];
				T5 = I[WS(is, 47)];
				T6 = T4 - T5;
				T338 = T4 + T5;
				T7 = I[WS(is, 15)];
				T8 = I[WS(is, 48)];
				T9 = T7 - T8;
				T339 = T7 + T8;
			}
			T273 = T271 - T272;
			T340 = T338 - T339;
			{
				E T10;
				E T146;
				E T385;
				E T386;
				T10 = KP707106781 * (T6 + T9);
				T11 = T3 - T10;
				T195 = T3 + T10;
				T146 = KP707106781 * (T6 - T9);
				T150 = T146 - T149;
				T228 = T149 + T146;
				T385 = T271 + T272;
				T386 = T338 + T339;
				T387 = T385 - T386;
				T431 = T385 + T386;
			}
		}
		{
			E T14;
			E T274;
			E T24;
			E T277;
			E T17;
			E T275;
			E T21;
			E T278;
			{
				E T12;
				E T13;
				E T22;
				E T23;
				T12 = I[WS(is, 8)];
				T13 = I[WS(is, 55)];
				T14 = T12 - T13;
				T274 = T12 + T13;
				T22 = I[WS(is, 7)];
				T23 = I[WS(is, 56)];
				T24 = T22 - T23;
				T277 = T22 + T23;
			}
			{
				E T15;
				E T16;
				E T19;
				E T20;
				T15 = I[WS(is, 40)];
				T16 = I[WS(is, 23)];
				T17 = T15 - T16;
				T275 = T15 + T16;
				T19 = I[WS(is, 24)];
				T20 = I[WS(is, 39)];
				T21 = T19 - T20;
				T278 = T19 + T20;
			}
			{
				E T18;
				E T25;
				E T416;
				E T417;
				T18 = FMA(KP382683432, T14, KP923879532 * T17);
				T25 = FNMS(KP382683432, T24, KP923879532 * T21);
				T26 = T18 - T25;
				T229 = T18 + T25;
				T416 = T274 + T275;
				T417 = T277 + T278;
				T418 = T416 - T417;
				T432 = T416 + T417;
			}
			{
				E T143;
				E T144;
				E T276;
				E T279;
				T143 = FNMS(KP382683432, T17, KP923879532 * T14);
				T144 = FMA(KP923879532, T24, KP382683432 * T21);
				T145 = T143 - T144;
				T196 = T143 + T144;
				T276 = T274 - T275;
				T279 = T277 - T278;
				T280 = KP707106781 * (T276 + T279);
				T337 = KP707106781 * (T276 - T279);
			}
		}
		{
			E T37;
			E T283;
			E T41;
			E T282;
			E T34;
			E T286;
			E T42;
			E T285;
			{
				E T35;
				E T36;
				E T39;
				E T40;
				T35 = I[WS(is, 36)];
				T36 = I[WS(is, 27)];
				T37 = T35 - T36;
				T283 = T35 + T36;
				T39 = I[WS(is, 4)];
				T40 = I[WS(is, 59)];
				T41 = T39 - T40;
				T282 = T39 + T40;
				{
					E T28;
					E T29;
					E T30;
					E T31;
					E T32;
					E T33;
					T28 = I[WS(is, 20)];
					T29 = I[WS(is, 43)];
					T30 = T28 - T29;
					T31 = I[WS(is, 11)];
					T32 = I[WS(is, 52)];
					T33 = T31 - T32;
					T34 = KP707106781 * (T30 - T33);
					T286 = T31 + T32;
					T42 = KP707106781 * (T30 + T33);
					T285 = T28 + T29;
				}
			}
			{
				E T38;
				E T43;
				E T388;
				E T389;
				T38 = T34 - T37;
				T43 = T41 - T42;
				T44 = FMA(KP831469612, T38, KP555570233 * T43);
				T140 = FNMS(KP555570233, T38, KP831469612 * T43);
				T388 = T282 + T283;
				T389 = T285 + T286;
				T390 = T388 - T389;
				T443 = T388 + T389;
			}
			{
				E T198;
				E T199;
				E T284;
				E T287;
				T198 = T41 + T42;
				T199 = T37 + T34;
				T200 = FNMS(KP195090322, T199, KP980785280 * T198);
				T225 = FMA(KP980785280, T199, KP195090322 * T198);
				T284 = T282 - T283;
				T287 = T285 - T286;
				T288 = FMA(KP382683432, T284, KP923879532 * T287);
				T334 = FNMS(KP382683432, T287, KP923879532 * T284);
			}
		}
		{
			E T54;
			E T293;
			E T58;
			E T292;
			E T51;
			E T290;
			E T59;
			E T289;
			{
				E T52;
				E T53;
				E T56;
				E T57;
				T52 = I[WS(is, 28)];
				T53 = I[WS(is, 35)];
				T54 = T52 - T53;
				T293 = T52 + T53;
				T56 = I[WS(is, 3)];
				T57 = I[WS(is, 60)];
				T58 = T56 - T57;
				T292 = T56 + T57;
				{
					E T45;
					E T46;
					E T47;
					E T48;
					E T49;
					E T50;
					T45 = I[WS(is, 12)];
					T46 = I[WS(is, 51)];
					T47 = T45 - T46;
					T48 = I[WS(is, 19)];
					T49 = I[WS(is, 44)];
					T50 = T48 - T49;
					T51 = KP707106781 * (T47 - T50);
					T290 = T48 + T49;
					T59 = KP707106781 * (T47 + T50);
					T289 = T45 + T46;
				}
			}
			{
				E T55;
				E T60;
				E T391;
				E T392;
				T55 = T51 - T54;
				T60 = T58 - T59;
				T61 = FNMS(KP555570233, T60, KP831469612 * T55);
				T141 = FMA(KP555570233, T55, KP831469612 * T60);
				T391 = T292 + T293;
				T392 = T289 + T290;
				T393 = T391 - T392;
				T444 = T391 + T392;
			}
			{
				E T201;
				E T202;
				E T291;
				E T294;
				T201 = T54 + T51;
				T202 = T58 + T59;
				T203 = FMA(KP195090322, T201, KP980785280 * T202);
				T226 = FNMS(KP195090322, T202, KP980785280 * T201);
				T291 = T289 - T290;
				T294 = T292 - T293;
				T295 = FNMS(KP382683432, T294, KP923879532 * T291);
				T335 = FMA(KP923879532, T294, KP382683432 * T291);
			}
		}
		{
			E T125;
			E T327;
			E T130;
			E T326;
			E T122;
			E T323;
			E T131;
			E T322;
			E T315;
			E T316;
			E T107;
			E T317;
			E T133;
			E T318;
			E T319;
			E T114;
			E T320;
			E T134;
			{
				E T123;
				E T124;
				E T128;
				E T129;
				T123 = I[WS(is, 30)];
				T124 = I[WS(is, 33)];
				T125 = T123 - T124;
				T327 = T123 + T124;
				T128 = I[WS(is, 1)];
				T129 = I[WS(is, 62)];
				T130 = T128 - T129;
				T326 = T128 + T129;
			}
			{
				E T116;
				E T117;
				E T118;
				E T119;
				E T120;
				E T121;
				T116 = I[WS(is, 14)];
				T117 = I[WS(is, 49)];
				T118 = T116 - T117;
				T119 = I[WS(is, 17)];
				T120 = I[WS(is, 46)];
				T121 = T119 - T120;
				T122 = KP707106781 * (T118 - T121);
				T323 = T119 + T120;
				T131 = KP707106781 * (T118 + T121);
				T322 = T116 + T117;
			}
			{
				E T103;
				E T106;
				E T110;
				E T113;
				{
					E T101;
					E T102;
					E T104;
					E T105;
					T101 = I[WS(is, 6)];
					T102 = I[WS(is, 57)];
					T103 = T101 - T102;
					T315 = T101 + T102;
					T104 = I[WS(is, 38)];
					T105 = I[WS(is, 25)];
					T106 = T104 - T105;
					T316 = T104 + T105;
				}
				T107 = FNMS(KP382683432, T106, KP923879532 * T103);
				T317 = T315 - T316;
				T133 = FMA(KP382683432, T103, KP923879532 * T106);
				{
					E T108;
					E T109;
					E T111;
					E T112;
					T108 = I[WS(is, 9)];
					T109 = I[WS(is, 54)];
					T110 = T108 - T109;
					T318 = T108 + T109;
					T111 = I[WS(is, 22)];
					T112 = I[WS(is, 41)];
					T113 = T111 - T112;
					T319 = T111 + T112;
				}
				T114 = FMA(KP923879532, T110, KP382683432 * T113);
				T320 = T318 - T319;
				T134 = FNMS(KP382683432, T110, KP923879532 * T113);
			}
			{
				E T115;
				E T126;
				E T403;
				E T404;
				T115 = T107 - T114;
				T126 = T122 - T125;
				T127 = T115 - T126;
				T173 = T126 + T115;
				T403 = T315 + T316;
				T404 = T318 + T319;
				T405 = T403 - T404;
				T438 = T403 + T404;
			}
			{
				E T406;
				E T407;
				E T132;
				E T135;
				T406 = T326 + T327;
				T407 = T322 + T323;
				T408 = T406 - T407;
				T437 = T406 + T407;
				T132 = T130 - T131;
				T135 = T133 - T134;
				T136 = T132 - T135;
				T174 = T132 + T135;
			}
			{
				E T213;
				E T214;
				E T321;
				E T324;
				T213 = T125 + T122;
				T214 = T133 + T134;
				T215 = T213 + T214;
				T249 = T214 - T213;
				T321 = KP707106781 * (T317 - T320);
				T324 = T322 - T323;
				T325 = T321 - T324;
				T363 = T324 + T321;
			}
			{
				E T328;
				E T329;
				E T216;
				E T217;
				T328 = T326 - T327;
				T329 = KP707106781 * (T317 + T320);
				T330 = T328 - T329;
				T364 = T328 + T329;
				T216 = T130 + T131;
				T217 = T107 + T114;
				T218 = T216 + T217;
				T250 = T216 - T217;
			}
		}
		{
			E T88;
			E T310;
			E T93;
			E T309;
			E T85;
			E T306;
			E T94;
			E T305;
			E T298;
			E T299;
			E T70;
			E T300;
			E T96;
			E T301;
			E T302;
			E T77;
			E T303;
			E T97;
			{
				E T86;
				E T87;
				E T91;
				E T92;
				T86 = I[WS(is, 34)];
				T87 = I[WS(is, 29)];
				T88 = T86 - T87;
				T310 = T86 + T87;
				T91 = I[WS(is, 2)];
				T92 = I[WS(is, 61)];
				T93 = T91 - T92;
				T309 = T91 + T92;
			}
			{
				E T79;
				E T80;
				E T81;
				E T82;
				E T83;
				E T84;
				T79 = I[WS(is, 18)];
				T80 = I[WS(is, 45)];
				T81 = T79 - T80;
				T82 = I[WS(is, 13)];
				T83 = I[WS(is, 50)];
				T84 = T82 - T83;
				T85 = KP707106781 * (T81 - T84);
				T306 = T82 + T83;
				T94 = KP707106781 * (T81 + T84);
				T305 = T79 + T80;
			}
			{
				E T66;
				E T69;
				E T73;
				E T76;
				{
					E T64;
					E T65;
					E T67;
					E T68;
					T64 = I[WS(is, 10)];
					T65 = I[WS(is, 53)];
					T66 = T64 - T65;
					T298 = T64 + T65;
					T67 = I[WS(is, 42)];
					T68 = I[WS(is, 21)];
					T69 = T67 - T68;
					T299 = T67 + T68;
				}
				T70 = FNMS(KP382683432, T69, KP923879532 * T66);
				T300 = T298 - T299;
				T96 = FMA(KP382683432, T66, KP923879532 * T69);
				{
					E T71;
					E T72;
					E T74;
					E T75;
					T71 = I[WS(is, 5)];
					T72 = I[WS(is, 58)];
					T73 = T71 - T72;
					T301 = T71 + T72;
					T74 = I[WS(is, 26)];
					T75 = I[WS(is, 37)];
					T76 = T74 - T75;
					T302 = T74 + T75;
				}
				T77 = FMA(KP923879532, T73, KP382683432 * T76);
				T303 = T301 - T302;
				T97 = FNMS(KP382683432, T73, KP923879532 * T76);
			}
			{
				E T78;
				E T89;
				E T396;
				E T397;
				T78 = T70 - T77;
				T89 = T85 - T88;
				T90 = T78 - T89;
				T170 = T89 + T78;
				T396 = T309 + T310;
				T397 = T305 + T306;
				T398 = T396 - T397;
				T434 = T396 + T397;
			}
			{
				E T399;
				E T400;
				E T95;
				E T98;
				T399 = T298 + T299;
				T400 = T301 + T302;
				T401 = T399 - T400;
				T435 = T399 + T400;
				T95 = T93 - T94;
				T98 = T96 - T97;
				T99 = T95 - T98;
				T171 = T95 + T98;
			}
			{
				E T206;
				E T207;
				E T304;
				E T307;
				T206 = T93 + T94;
				T207 = T70 + T77;
				T208 = T206 + T207;
				T246 = T206 - T207;
				T304 = KP707106781 * (T300 - T303);
				T307 = T305 - T306;
				T308 = T304 - T307;
				T360 = T307 + T304;
			}
			{
				E T311;
				E T312;
				E T209;
				E T210;
				T311 = T309 - T310;
				T312 = KP707106781 * (T300 + T303);
				T313 = T311 - T312;
				T361 = T311 + T312;
				T209 = T88 + T85;
				T210 = T96 + T97;
				T211 = T209 + T210;
				T247 = T210 - T209;
			}
		}
		{
			E T451;
			E T455;
			E T454;
			E T456;
			{
				E T449;
				E T450;
				E T452;
				E T453;
				T449 = T431 + T432;
				T450 = T443 + T444;
				T451 = T449 - T450;
				T455 = T449 + T450;
				T452 = T434 + T435;
				T453 = T437 + T438;
				T454 = T452 - T453;
				T456 = T452 + T453;
			}
			O[WS(os, 16)] = FNMS(KP765366864, T454, KP1_847759065 * T451);
			O[0] = KP2_000000000 * (T455 + T456);
			O[WS(os, 48)] = FMA(KP765366864, T451, KP1_847759065 * T454);
			O[WS(os, 32)] = KP1_414213562 * (T455 - T456);
		}
		{
			E T433;
			E T445;
			E T440;
			E T442;
			E T436;
			E T439;
			T433 = T431 - T432;
			T445 = T443 - T444;
			T436 = T434 - T435;
			T439 = T437 - T438;
			T440 = KP707106781 * (T436 + T439);
			T442 = KP707106781 * (T436 - T439);
			{
				E T441;
				E T446;
				E T447;
				E T448;
				T441 = T433 - T440;
				T446 = T442 - T445;
				O[WS(os, 24)] = FNMS(KP1_111140466, T446, KP1_662939224 * T441);
				O[WS(os, 40)] = FMA(KP1_662939224, T446, KP1_111140466 * T441);
				T447 = T433 + T440;
				T448 = T445 + T442;
				O[WS(os, 8)] = FNMS(KP390180644, T448, KP1_961570560 * T447);
				O[WS(os, 56)] = FMA(KP1_961570560, T448, KP390180644 * T447);
			}
		}
		{
			E T395;
			E T423;
			E T419;
			E T426;
			E T410;
			E T427;
			E T414;
			E T424;
			E T394;
			E T415;
			T394 = KP707106781 * (T390 + T393);
			T395 = T387 - T394;
			T423 = T387 + T394;
			T415 = KP707106781 * (T390 - T393);
			T419 = T415 - T418;
			T426 = T418 + T415;
			{
				E T402;
				E T409;
				E T412;
				E T413;
				T402 = FMA(KP382683432, T398, KP923879532 * T401);
				T409 = FNMS(KP382683432, T408, KP923879532 * T405);
				T410 = T402 - T409;
				T427 = T402 + T409;
				T412 = FNMS(KP382683432, T401, KP923879532 * T398);
				T413 = FMA(KP923879532, T408, KP382683432 * T405);
				T414 = T412 - T413;
				T424 = T412 + T413;
			}
			{
				E T411;
				E T420;
				E T429;
				E T430;
				T411 = T395 - T410;
				T420 = T414 - T419;
				O[WS(os, 20)] = FNMS(KP942793473, T420, KP1_763842528 * T411);
				O[WS(os, 44)] = FMA(KP1_763842528, T420, KP942793473 * T411);
				T429 = T423 - T424;
				T430 = T427 - T426;
				O[WS(os, 28)] = FNMS(KP1_268786568, T430, KP1_546020906 * T429);
				O[WS(os, 36)] = FMA(KP1_268786568, T429, KP1_546020906 * T430);
			}
			{
				E T421;
				E T422;
				E T425;
				E T428;
				T421 = T395 + T410;
				T422 = T419 + T414;
				O[WS(os, 12)] = FNMS(KP580569354, T422, KP1_913880671 * T421);
				O[WS(os, 52)] = FMA(KP1_913880671, T422, KP580569354 * T421);
				T425 = T423 + T424;
				T428 = T426 + T427;
				O[WS(os, 4)] = FNMS(KP196034280, T428, KP1_990369453 * T425);
				O[WS(os, 60)] = FMA(KP196034280, T425, KP1_990369453 * T428);
			}
		}
		{
			E T359;
			E T377;
			E T373;
			E T378;
			E T366;
			E T380;
			E T370;
			E T381;
			{
				E T357;
				E T358;
				E T371;
				E T372;
				T357 = T273 + T280;
				T358 = T334 + T335;
				T359 = T357 - T358;
				T377 = T357 + T358;
				T371 = FNMS(KP195090322, T360, KP980785280 * T361);
				T372 = FMA(KP195090322, T363, KP980785280 * T364);
				T373 = T371 - T372;
				T378 = T371 + T372;
			}
			{
				E T362;
				E T365;
				E T368;
				E T369;
				T362 = FMA(KP980785280, T360, KP195090322 * T361);
				T365 = FNMS(KP195090322, T364, KP980785280 * T363);
				T366 = T362 - T365;
				T380 = T362 + T365;
				T368 = T288 + T295;
				T369 = T340 + T337;
				T370 = T368 - T369;
				T381 = T369 + T368;
			}
			{
				E T367;
				E T374;
				E T383;
				E T384;
				T367 = T359 + T366;
				T374 = T370 + T373;
				O[WS(os, 14)] = FNMS(KP673779706, T374, KP1_883088130 * T367);
				O[WS(os, 50)] = FMA(KP673779706, T367, KP1_883088130 * T374);
				T383 = T377 + T378;
				T384 = T381 + T380;
				O[WS(os, 2)] = FNMS(KP098135348, T384, KP1_997590912 * T383);
				O[WS(os, 62)] = FMA(KP1_997590912, T384, KP098135348 * T383);
			}
			{
				E T375;
				E T376;
				E T379;
				E T382;
				T375 = T359 - T366;
				T376 = T373 - T370;
				O[WS(os, 18)] = FNMS(KP855110186, T376, KP1_807978586 * T375);
				O[WS(os, 46)] = FMA(KP855110186, T375, KP1_807978586 * T376);
				T379 = T377 - T378;
				T382 = T380 - T381;
				O[WS(os, 30)] = FNMS(KP1_343117909, T382, KP1_481902250 * T379);
				O[WS(os, 34)] = FMA(KP1_481902250, T382, KP1_343117909 * T379);
			}
		}
		{
			E T63;
			E T159;
			E T155;
			E T160;
			E T138;
			E T162;
			E T152;
			E T163;
			{
				E T27;
				E T62;
				E T153;
				E T154;
				T27 = T11 - T26;
				T62 = T44 - T61;
				T63 = T27 - T62;
				T159 = T27 + T62;
				T153 = FNMS(KP471396736, T90, KP881921264 * T99);
				T154 = FMA(KP471396736, T127, KP881921264 * T136);
				T155 = T153 - T154;
				T160 = T153 + T154;
			}
			{
				E T100;
				E T137;
				E T142;
				E T151;
				T100 = FMA(KP881921264, T90, KP471396736 * T99);
				T137 = FNMS(KP471396736, T136, KP881921264 * T127);
				T138 = T100 - T137;
				T162 = T100 + T137;
				T142 = T140 - T141;
				T151 = T145 - T150;
				T152 = T142 - T151;
				T163 = T151 + T142;
			}
			{
				E T139;
				E T156;
				E T165;
				E T166;
				T139 = T63 + T138;
				T156 = T152 + T155;
				O[WS(os, 11)] = FNMS(KP533425514, T156, KP1_927552131 * T139);
				O[WS(os, 53)] = FMA(KP533425514, T139, KP1_927552131 * T156);
				T165 = T159 + T160;
				T166 = T163 + T162;
				O[WS(os, 5)] = FNMS(KP244821350, T166, KP1_984959069 * T165);
				O[WS(os, 59)] = FMA(KP1_984959069, T166, KP244821350 * T165);
			}
			{
				E T157;
				E T158;
				E T161;
				E T164;
				T157 = T63 - T138;
				T158 = T155 - T152;
				O[WS(os, 21)] = FNMS(KP985796384, T158, KP1_740173982 * T157);
				O[WS(os, 43)] = FMA(KP985796384, T157, KP1_740173982 * T158);
				T161 = T159 - T160;
				T164 = T162 - T163;
				O[WS(os, 27)] = FNMS(KP1_230463181, T164, KP1_576692855 * T161);
				O[WS(os, 37)] = FMA(KP1_576692855, T164, KP1_230463181 * T161);
			}
		}
		{
			E T205;
			E T235;
			E T231;
			E T238;
			E T220;
			E T239;
			E T224;
			E T236;
			{
				E T197;
				E T204;
				E T227;
				E T230;
				T197 = T195 + T196;
				T204 = T200 + T203;
				T205 = T197 - T204;
				T235 = T197 + T204;
				T227 = T225 + T226;
				T230 = T228 + T229;
				T231 = T227 - T230;
				T238 = T230 + T227;
			}
			{
				E T212;
				E T219;
				E T222;
				E T223;
				T212 = FMA(KP098017140, T208, KP995184726 * T211);
				T219 = FNMS(KP098017140, T218, KP995184726 * T215);
				T220 = T212 - T219;
				T239 = T212 + T219;
				T222 = FNMS(KP098017140, T211, KP995184726 * T208);
				T223 = FMA(KP995184726, T218, KP098017140 * T215);
				T224 = T222 - T223;
				T236 = T222 + T223;
			}
			{
				E T221;
				E T232;
				E T241;
				E T242;
				T221 = T205 - T220;
				T232 = T224 - T231;
				O[WS(os, 17)] = FNMS(KP810482628, T232, KP1_828419511 * T221);
				O[WS(os, 47)] = FMA(KP1_828419511, T232, KP810482628 * T221);
				T241 = T235 - T236;
				T242 = T239 - T238;
				O[WS(os, 31)] = FNMS(KP1_379081089, T242, KP1_448494165 * T241);
				O[WS(os, 33)] = FMA(KP1_379081089, T241, KP1_448494165 * T242);
			}
			{
				E T233;
				E T234;
				E T237;
				E T240;
				T233 = T205 + T220;
				T234 = T231 + T224;
				O[WS(os, 15)] = FNMS(KP719790073, T234, KP1_865985597 * T233);
				O[WS(os, 49)] = FMA(KP1_865985597, T234, KP719790073 * T233);
				T237 = T235 + T236;
				T240 = T238 + T239;
				O[WS(os, 1)] = FNMS(KP049082457, T240, KP1_999397637 * T237);
				O[WS(os, 63)] = FMA(KP049082457, T237, KP1_999397637 * T240);
			}
		}
		{
			E T297;
			E T349;
			E T345;
			E T350;
			E T332;
			E T352;
			E T342;
			E T353;
			{
				E T281;
				E T296;
				E T343;
				E T344;
				T281 = T273 - T280;
				T296 = T288 - T295;
				T297 = T281 - T296;
				T349 = T281 + T296;
				T343 = FNMS(KP555570233, T308, KP831469612 * T313);
				T344 = FMA(KP555570233, T325, KP831469612 * T330);
				T345 = T343 - T344;
				T350 = T343 + T344;
			}
			{
				E T314;
				E T331;
				E T336;
				E T341;
				T314 = FMA(KP831469612, T308, KP555570233 * T313);
				T331 = FNMS(KP555570233, T330, KP831469612 * T325);
				T332 = T314 - T331;
				T352 = T314 + T331;
				T336 = T334 - T335;
				T341 = T337 - T340;
				T342 = T336 - T341;
				T353 = T341 + T336;
			}
			{
				E T333;
				E T346;
				E T355;
				E T356;
				T333 = T297 + T332;
				T346 = T342 + T345;
				O[WS(os, 10)] = FNMS(KP485960359, T346, KP1_940062506 * T333);
				O[WS(os, 54)] = FMA(KP485960359, T333, KP1_940062506 * T346);
				T355 = T349 + T350;
				T356 = T353 + T352;
				O[WS(os, 6)] = FNMS(KP293460948, T356, KP1_978353019 * T355);
				O[WS(os, 58)] = FMA(KP1_978353019, T356, KP293460948 * T355);
			}
			{
				E T347;
				E T348;
				E T351;
				E T354;
				T347 = T297 - T332;
				T348 = T345 - T342;
				O[WS(os, 22)] = FNMS(KP1_028205488, T348, KP1_715457220 * T347);
				O[WS(os, 42)] = FMA(KP1_028205488, T347, KP1_715457220 * T348);
				T351 = T349 - T350;
				T354 = T352 - T353;
				O[WS(os, 26)] = FNMS(KP1_191398608, T354, KP1_606415062 * T351);
				O[WS(os, 38)] = FMA(KP1_606415062, T354, KP1_191398608 * T351);
			}
		}
		{
			E T169;
			E T187;
			E T183;
			E T188;
			E T176;
			E T190;
			E T180;
			E T191;
			{
				E T167;
				E T168;
				E T181;
				E T182;
				T167 = T11 + T26;
				T168 = T140 + T141;
				T169 = T167 - T168;
				T187 = T167 + T168;
				T181 = FNMS(KP290284677, T170, KP956940335 * T171);
				T182 = FMA(KP290284677, T173, KP956940335 * T174);
				T183 = T181 - T182;
				T188 = T181 + T182;
			}
			{
				E T172;
				E T175;
				E T178;
				E T179;
				T172 = FMA(KP956940335, T170, KP290284677 * T171);
				T175 = FNMS(KP290284677, T174, KP956940335 * T173);
				T176 = T172 - T175;
				T190 = T172 + T175;
				T178 = T44 + T61;
				T179 = T150 + T145;
				T180 = T178 - T179;
				T191 = T179 + T178;
			}
			{
				E T177;
				E T184;
				E T193;
				E T194;
				T177 = T169 + T176;
				T184 = T180 + T183;
				O[WS(os, 13)] = FNMS(KP627363480, T184, KP1_899056361 * T177);
				O[WS(os, 51)] = FMA(KP627363480, T177, KP1_899056361 * T184);
				T193 = T187 + T188;
				T194 = T191 + T190;
				O[WS(os, 3)] = FNMS(KP147129127, T194, KP1_994580913 * T193);
				O[WS(os, 61)] = FMA(KP1_994580913, T194, KP147129127 * T193);
			}
			{
				E T185;
				E T186;
				E T189;
				E T192;
				T185 = T169 - T176;
				T186 = T183 - T180;
				O[WS(os, 19)] = FNMS(KP899222659, T186, KP1_786448602 * T185);
				O[WS(os, 45)] = FMA(KP899222659, T185, KP1_786448602 * T186);
				T189 = T187 - T188;
				T192 = T190 - T191;
				O[WS(os, 29)] = FNMS(KP1_306345685, T192, KP1_514417693 * T189);
				O[WS(os, 35)] = FMA(KP1_514417693, T192, KP1_306345685 * T189);
			}
		}
		{
			E T245;
			E T263;
			E T259;
			E T266;
			E T252;
			E T267;
			E T256;
			E T264;
			{
				E T243;
				E T244;
				E T257;
				E T258;
				T243 = T195 - T196;
				T244 = T225 - T226;
				T245 = T243 - T244;
				T263 = T243 + T244;
				T257 = T200 - T203;
				T258 = T229 - T228;
				T259 = T257 - T258;
				T266 = T258 + T257;
			}
			{
				E T248;
				E T251;
				E T254;
				E T255;
				T248 = FMA(KP634393284, T246, KP773010453 * T247);
				T251 = FNMS(KP634393284, T250, KP773010453 * T249);
				T252 = T248 - T251;
				T267 = T248 + T251;
				T254 = FNMS(KP634393284, T247, KP773010453 * T246);
				T255 = FMA(KP773010453, T250, KP634393284 * T249);
				T256 = T254 - T255;
				T264 = T254 + T255;
			}
			{
				E T253;
				E T260;
				E T269;
				E T270;
				T253 = T245 - T252;
				T260 = T256 - T259;
				O[WS(os, 23)] = FNMS(KP1_069995239, T260, KP1_689707130 * T253);
				O[WS(os, 41)] = FMA(KP1_689707130, T260, KP1_069995239 * T253);
				T269 = T263 - T264;
				T270 = T267 - T266;
				O[WS(os, 25)] = FNMS(KP1_151616382, T270, KP1_635169626 * T269);
				O[WS(os, 39)] = FMA(KP1_151616382, T269, KP1_635169626 * T270);
			}
			{
				E T261;
				E T262;
				E T265;
				E T268;
				T261 = T245 + T252;
				T262 = T259 + T256;
				O[WS(os, 9)] = FNMS(KP438202480, T262, KP1_951404260 * T261);
				O[WS(os, 55)] = FMA(KP1_951404260, T262, KP438202480 * T261);
				T265 = T263 + T264;
				T268 = T266 + T267;
				O[WS(os, 7)] = FNMS(KP341923777, T268, KP1_970555284 * T265);
				O[WS(os, 57)] = FMA(KP341923777, T265, KP1_970555284 * T268);
			}
		}
	}
}


template void dct_type2_64<>(const double*, double*, int, int, int, int, int);
template void dct_type2_64<>(const float*, float*, int, int, int, int, int);
template void dct_type2_64<>(const double*, double*, size_t, size_t, size_t, size_t, size_t);
template void dct_type2_64<>(const float*, float*, size_t, size_t, size_t, size_t, size_t);
