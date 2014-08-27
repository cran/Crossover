package org.mutoss.gui;

import java.awt.Font;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.io.File;
import java.io.IOException;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;

import org.af.commons.io.FileTools;
import org.mutoss.config.ClassConfig;
import org.mutoss.config.Configuration;

import com.jgoodies.forms.layout.CellConstraints;
import com.jgoodies.forms.layout.FormLayout;

public class DesignInputPanel extends JPanel implements KeyListener, ActionListener {
	
	JTextField jtTitle = new JTextField();
	JTextField jtReference = new JTextField();
	JButton ok = new JButton("Ready");
	JButton loadFile = new JButton("Load File");
	JButton loadRObject = new JButton("Load from R");
	JButton save = new JButton("Save to my archive of designs");
	JTextArea jta;
	JLabel label = new JLabel();
	JComboBox jcbRows = new JComboBox(new String[] {"periods", "sequences"});
	CrossoverGUI gui;
	
	ClassConfig ac = new ClassConfig(Configuration.getInstance(), DesignInputPanel.class);
	
	public DesignInputPanel(CrossoverGUI gui) {
		this.gui = gui;
		String cols = "5dlu, fill:min:grow, 5dlu, fill:min:grow, 5dlu,";
        String rows = "5dlu, pref, 5dlu, fill:pref:grow, 5dlu, pref, 5dlu";
        
        FormLayout layout = new FormLayout(cols, rows);
        layout.setColumnGroups(new int[][]{ {2, 4} });

        setLayout(layout);
        CellConstraints cc = new CellConstraints();
		
		int row = 2;
    	
    	add(new JLabel("Please enter your design here:"), cc.xy(2, row));
        add(new JLabel(""), cc.xy(4, row));
		
        row+=2;
        
		jta = new JTextArea("");
		//jta.setBorder(LineBorder.createBlackLineBorder());
		jta.setFont(new Font("Monospaced", Font.PLAIN, 12));
		jta.setLineWrap(false);		
		jta.setMargin(new Insets(4,4,4,4));
		jta.addKeyListener(this);
		
		add(new JScrollPane(jta), cc.xy(2, row));
		add(new JScrollPane(getRightSidePanel()), cc.xy(4, row));
		
		row+=2;
		
		loadFile.addActionListener(this);
		add(loadFile, cc.xy(2, row));
		save.addActionListener(this);
		add(save, cc.xy(4, row));
		
		loadDefaults();
		
	}
	
	public JPanel getRightSidePanel() {
		JPanel panel = new JPanel();
		String cols = "5dlu, pref, 5dlu, fill:min:grow, 5dlu";
        String rows = "5dlu, pref, 5dlu, pref, 5dlu, pref, 5dlu, pref, 5dlu";
        
        panel.setLayout(new FormLayout(cols, rows));
        CellConstraints cc = new CellConstraints();

        int row = 2;

        panel.add(new JLabel("Title"), cc.xy(2, row));
        panel.add(jtTitle, cc.xy(4, row));
        
        row+=2;

        panel.add(new JLabel("Reference"), cc.xy(2, row));
        panel.add(jtReference, cc.xy(4, row));

        row+=2;

        jcbRows.addActionListener(this);
        panel.add(new JLabel("Rows represent"), cc.xy(2, row));
        panel.add(jcbRows, cc.xy(4, row));
        
        row+=2;

        panel.add(label, cc.xyw(2, row, 3));
        
        return panel;
	}
	
	public void loadDefaults() {
		jcbRows.setSelectedIndex(ac.getIntProperty("CVPattern", 0));
	}
    
	public void saveDefaults() {
		ac.setIntProperty("CVPattern", jcbRows.getSelectedIndex());
	}
	
	public boolean transpose() {
		return !jcbRows.getSelectedItem().equals("periods");
	}
	
	public void addActionListener(ActionListener al) {
		ok.addActionListener(al);
	}
	
	public Design getDesign() {
		return null;
	}

	public void keyPressed(KeyEvent e) {}

	public void keyReleased(KeyEvent e) {
		checkDesign();
	}

	public void keyTyped(KeyEvent e) {	}

	public void actionPerformed(ActionEvent e) {
		saveDefaults();
		if (e.getSource() == loadFile) {
			JFileChooser fc = new JFileChooser(Configuration.getInstance().getClassProperty(this.getClass(), "DesignDirectory"));		
			fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
	        int returnVal = fc.showOpenDialog(this);
	        if (returnVal == JFileChooser.APPROVE_OPTION) {
	        	File f = fc.getSelectedFile();
	        	Configuration.getInstance().setClassProperty(this.getClass(), "DesignDirectory", f.getParent());
	        	try {
					String text = FileTools.readFileAsString(f);
					jta.setText(text);
				} catch (IOException e1) {
					JOptionPane.showMessageDialog(this, "File could not be opened:\n"+e1.getMessage(), "Error opening file", JOptionPane.ERROR_MESSAGE);
				}	        	
	        }
		} else if (e.getSource() == jcbRows) {
			checkDesign();
		} else if (e.getSource() == save) {
			try {
				saveDesign();
				gui.dac.addEnteredDesign(new Design(jtTitle.getText(), ".newDesign", jtReference.getText()));
			} catch (Exception error) {
				JOptionPane.showMessageDialog(this, "Design not valid.\nPlease correct.", "Design not valid", JOptionPane.ERROR_MESSAGE);
			}
		}
	}
	
	/**
	 * Saves the design as .newDesign or throws error.
	 * @return
	 * @throws Exception 
	 */
	private void saveDesign() throws Exception {
		String input = jta.getText();
		RControl.getR().evalVoid(".con <- textConnection(\""+input+"\")");
		RControl.getR().eval(".newDesign <- try(as.matrix(read.table(.con, header = FALSE)), silent=TRUE)");
		if (transpose()) {
			RControl.getR().eval(".newDesign <- t(.newDesign)");
		}
		if (RControl.getR().eval("(\"try-error\" %in% class(.newDesign))").asRLogical().getData()[0]) throw new Exception();
	}

	private void checkDesign() {
		try {
			saveDesign();
			int[] dim = RControl.getR().eval("dim(.newDesign)").asRInteger().getData(); 
			int t = RControl.getR().eval("length(levels(as.factor(.newDesign)))").asRInteger().getData()[0];	
			label.setText("p = "+dim[0]+", s = "+dim[1]+", t = "+t);
		} catch (Exception error) {
			label.setText("Design is not valid, please correct it.");
		}
	}
	
}
