/**
 * texProseToHtml(raw)
 *
 * Pre-processes LaTeX text-mode commands in prose strings so that they render
 * correctly in HTML, before KaTeX handles any embedded math delimiters.
 *
 * Handled commands (outside math mode):
 *   \texttt{X}        → <code>X</code>
 *   \textbf{X}        → <strong>X</strong>
 *   \textit{X}        → <em>X</em>
 *   \emph{X}          → <em>X</em>
 *   \href{URL}{TEXT}  → <a href="URL">TEXT</a>  (http(s) URLs open in new tab)
 *   \url{URL}         → <a href="URL">URL</a>
 *
 * Math spans ($...$, \(...\), \[...\], $$...$$) are skipped so their
 * interior is not altered.
 */

// Match a balanced brace group: {…}  — handles one level of nesting.
// For the prose content in this catalog, one level is sufficient.
function extractBraceGroup(s: string, start: number): { content: string; end: number } | null {
  if (s[start] !== '{') return null;
  let depth = 0;
  let i = start;
  while (i < s.length) {
    if (s[i] === '{') depth++;
    else if (s[i] === '}') {
      depth--;
      if (depth === 0) return { content: s.slice(start + 1, i), end: i + 1 };
    }
    i++;
  }
  return null; // unbalanced
}

function isHttpUrl(url: string): boolean {
  return /^https?:\/\//i.test(url.trim());
}

export function texProseToHtml(raw: string): string {
  // We walk through the string character by character, skipping math spans so
  // we don't corrupt their content.  Outside math spans we replace text-mode
  // LaTeX commands with HTML equivalents.

  const mathOpeners: Array<{ open: string; close: string }> = [
    { open: '$$', close: '$$' },
    { open: '$', close: '$' },
    { open: '\\[', close: '\\]' },
    { open: '\\(', close: '\\)' },
  ];

  let out = '';
  let i = 0;
  const len = raw.length;

  while (i < len) {
    // Check if we're at the start of a math span.
    let mathMatch: { open: string; close: string } | null = null;
    for (const m of mathOpeners) {
      if (raw.startsWith(m.open, i)) {
        mathMatch = m;
        break;
      }
    }

    if (mathMatch) {
      // Emit the math span verbatim (KaTeX handles it later).
      const closeIdx = raw.indexOf(mathMatch.close, i + mathMatch.open.length);
      if (closeIdx === -1) {
        // No closing delimiter found — emit rest verbatim.
        out += raw.slice(i);
        break;
      }
      out += raw.slice(i, closeIdx + mathMatch.close.length);
      i = closeIdx + mathMatch.close.length;
      continue;
    }

    // Check for \texttt, \textbf, \textit, \emph, \href, \url
    if (raw[i] === '\\') {
      let handled = false;

      const commandPatterns: Array<{ name: string; tag: string; attrs?: string }> = [
        { name: 'texttt', tag: 'code' },
        { name: 'textbf', tag: 'strong' },
        { name: 'textit', tag: 'em' },
        { name: 'emph',   tag: 'em' },
      ];

      for (const { name, tag, attrs } of commandPatterns) {
        if (raw.startsWith('\\' + name + '{', i) || raw.startsWith('\\' + name + ' {', i)) {
          // Find the opening brace.
          let braceStart = i + 1 + name.length;
          while (braceStart < len && raw[braceStart] === ' ') braceStart++;
          const group = extractBraceGroup(raw, braceStart);
          if (group) {
            // Recursively apply transform to the group content.
            const inner = texProseToHtml(group.content);
            const attrStr = attrs ? ' ' + attrs : '';
            out += `<${tag}${attrStr}>${inner}</${tag}>`;
            i = group.end;
            handled = true;
            break;
          }
        }
      }

      if (!handled && raw.startsWith('\\href{', i)) {
        // \href{URL}{TEXT}
        const urlGroup = extractBraceGroup(raw, i + 5);
        if (urlGroup) {
          const textGroup = extractBraceGroup(raw, urlGroup.end);
          if (textGroup) {
            const url = urlGroup.content.trim();
            const text = texProseToHtml(textGroup.content);
            const targetAttr = isHttpUrl(url) ? ' target="_blank" rel="noopener noreferrer"' : '';
            out += `<a href="${url}"${targetAttr}>${text}</a>`;
            i = textGroup.end;
            handled = true;
          }
        }
      }

      if (!handled && raw.startsWith('\\url{', i)) {
        // \url{URL}
        const urlGroup = extractBraceGroup(raw, i + 4);
        if (urlGroup) {
          const url = urlGroup.content.trim();
          const targetAttr = isHttpUrl(url) ? ' target="_blank" rel="noopener noreferrer"' : '';
          out += `<a href="${url}"${targetAttr}>${url}</a>`;
          i = urlGroup.end;
          handled = true;
        }
      }

      if (!handled) {
        out += raw[i];
        i++;
      }
      continue;
    }

    out += raw[i];
    i++;
  }

  return out;
}
