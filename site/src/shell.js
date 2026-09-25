/* Site shell: tab routing (bare #hash tokens) and the export panel that replaces file downloads. */
(function () {
  'use strict';
  const TABS = ['explore', 'analysis', 'compare', 'about'];
  const listeners = [];
  window.Site = {
    current: 'explore',
    onTab(fn) { listeners.push(fn); },
    show(tab) {
      if (!TABS.includes(tab)) tab = 'explore';
      TABS.forEach(t => { document.getElementById(t).hidden = t !== tab; });
      document.querySelectorAll('[data-tab]').forEach(a => a.classList.toggle('on', a.dataset.tab === tab));
      this.current = tab;
      if (location.hash.slice(1) !== tab) history.replaceState(null, '', '#' + tab);
      listeners.forEach(fn => fn(tab));
    },
    // kind: 'image' (src = data URL) or 'text'
    exportPanel({ kind, src, text, name }) {
      const box = document.getElementById('xport'), body = document.getElementById('xport-body');
      document.getElementById('xport-title').textContent = name;
      body.innerHTML = '';
      if (kind === 'image') {
        document.getElementById('xport-hint').textContent = 'Right-click (or long-press) the image and choose “Save image as…”.';
        const img = new Image(); img.src = src; img.alt = name; body.appendChild(img);
      } else {
        document.getElementById('xport-hint').textContent = 'Copy the text below and save it as ' + name + '.';
        const ta = document.createElement('textarea'); ta.value = text; ta.readOnly = true; ta.id = 'xport-text'; body.appendChild(ta);
        const b = document.createElement('button'); b.className = 'btn'; b.textContent = 'Copy';
        b.onclick = () => navigator.clipboard.writeText(text).then(() => { b.textContent = 'Copied'; }, () => { ta.select(); });
        body.appendChild(b);
      }
      // real download button when not inside a sandboxed frame (e.g. on GitHub Pages)
      let top = false; try { top = window.self === window.top; } catch (e) { top = false; }
      if (top) {
        const a = document.createElement('a'); a.className = 'btn'; a.textContent = 'Download ' + name; a.download = name;
        a.href = kind === 'image' ? src : URL.createObjectURL(new Blob([text], { type: 'text/plain' }));
        a.style.display = 'inline-block'; a.style.marginTop = '8px'; a.style.textDecoration = 'none';
        body.appendChild(a);
        document.getElementById('xport-hint').textContent = 'Download the file, or ' + (kind === 'image' ? 'right-click the image to save it.' : 'copy the text below.');
      }
      box.hidden = false;
    }
  };
  document.addEventListener('click', e => {
    const a = e.target.closest('[data-tab]');
    if (a) { e.preventDefault(); Site.show(a.dataset.tab); window.scrollTo(0, 0); }
  });
  document.getElementById('xport-close').onclick = () => { document.getElementById('xport').hidden = true; };
  document.getElementById('xport').addEventListener('click', e => { if (e.target.id === 'xport') e.currentTarget.hidden = true; });
  addEventListener('keydown', e => { if (e.key === 'Escape') document.getElementById('xport').hidden = true; });
  addEventListener('hashchange', () => Site.show(location.hash.slice(1)));
  addEventListener('DOMContentLoaded', () => Site.show(location.hash.slice(1) || 'explore'));
})();
