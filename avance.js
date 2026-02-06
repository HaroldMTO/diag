function avance(val)
{
	var i,j,k,step0,steps,imgs,titles,maps,s;

	maps = document.getElementsByName("map");
	imgs = document.getElementsByName("fig");
	titles = document.getElementsByName("title");
	steps = document.getElementsByName("step");

	if (val == 0) {
		steps[0].options[0].selected = true;
		steps[0].options[0].selectedIndex = 0;
		return;
	}

	step0 = steps[0];
	i = step0.selectedIndex;

	s = step0.options[i].text;

	switch (val) {
	case -2:
		j = 0;
		break;
	case 2:
		j = step0.length-1;
		break;
	default:
		if (i+val < 0) j = 0;
		else if (i+val >= step0.length) j = step0.length-1;
		else j = i+val;
		break;
	}

	if (j == i) return;

	for (k=0;k<imgs.length;k++) imgs[k].src = maps[k].options[j].text;

	for (k=0;k<titles.length;k++) titles[k].innerHTML = step0.options[j].text;

	steps[0].options[i].selected = false;
	steps[0].options[j].selected = true;
	steps[0].selectedIndex = j;
}

function avanced(val)
{
	var i,j,k,l,dom0,doms,imgs,maps,s;

	maps = document.getElementsByName("map");
	imgs = document.getElementsByName("fig");
	doms = document.getElementsByName("dom");

	if (val == 0) {
		doms[0].options[0].selected = true;
		doms[0].options[0].selectedIndex = 0;
		return;
	}

	dom0 = doms[0];
	i = dom0.selectedIndex;

	s = dom0.options[i].text;

	if (i+val < 0) {
		j = 0;
	} else if (i+val >= dom0.length) {
		j = dom0.length-1;
	} else {
		j = i+val;
	}

	if (j == i) return;

	for (k=0;k<imgs.length;k++) imgs[k].src = imgs[k].src.replace(s,dom0.options[j].text);

	for (k=0;k<maps.length;k++) {
		for (l=0;l<maps[k].length;l++) {
			maps[k].options[l].text = maps[k].options[l].text.replace(s,dom0.options[j].text);
		}
	}

	doms[0].options[i].selected = false;
	doms[0].options[j].selected = true;
	doms[0].selectedIndex = j;
}
